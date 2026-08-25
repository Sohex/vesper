# The ecosystem tier read at source: ECOGEM against MARBL

WORLDBUILDING CONTEXT, stated first because this document borrows vocabulary
from several real disciplines: **Vesper is a fictional planet and this is
engineering work on the simulation of it.** Plankton, chlorophyll, photic zone,
productivity, Redfield and diatom name modelled quantities and modelled
parameters of an invented world throughout. Nothing here is Earth science and
nothing here is biology.

Read 2026-08-22 from `vendor/cgenie/genie-ecogem` and `references/marbl`, at
line level, without compiling or running either. This is the scoping note OCN-4
asks for: what each candidate for the ecosystem tier is, what it would cost to
drive from an offline circulation, and which of the two answers to which of the
row's questions. Every number below is a source constant or a count of source
lines, measured on the trees as vendored on this date.

The row's own summary was taken as a hypothesis and checked, not as a premise.
Most of it holds. Six parts of it need correcting, and one of the six inverts a
verdict.

---

## 1. Provenance, licence, version

**ECOGEM** is not a separate distribution. It is `genie-ecogem`, one of twenty
`genie-*` directories inside cGENIE.muffin, vendored at `vendor/cgenie/` as a
git subtree. MIT, `vendor/cgenie/LICENSE:1`, "Copyright (c) 2020 Andrew
Ridgwell". The tree's own version marker is
`vendor/cgenie/tagged_releases.readme:5`, "LATEST == v0.9.50", and the tagged
release for ECOGEM's own publication is at `tagged_releases.readme`, 0.9.0 and
0.9.1 against Ward et al. (2018) GMD. Seven Fortran files, 4395 lines:

    ecogem.f90              1011   the loop, and the BIOGEM coupling arrays
    ecogem_data_netCDF.f90  1036   output and the netCDF restart
    ecogem_data.f90          915   namelist load, trait construction, restart read
    ecogem_box.f90           583   uptake, photosynthesis, grazing, limitation
    ecogem_lib.f90           391   module state and the whole of ini_ecogem_nml
    initialise_ecogem.f90    398   index construction and allocation
    end_ecogem.f90            61   shutdown

**MARBL** is a standalone library, `references/marbl/`, BSD-style from UCAR,
`references/marbl/LICENSE.txt:1`. `references/INDEX.md:57` records the
provenance: `marbl-ecosys/MARBL` commit `403f6c50bba305d59910031dae33ff0629fb8602`.
The documentation version string is `docs/src/conf.py:57`, `cesm2.1`, with the
theme offering `latest` and `cesm2.1` side by side at `conf.py:99`. 35 Fortran
files under `src/`, plus `MARBL_tools/` (14 Python modules), `defaults/`
(7 YAML settings files and a JSON form, 11943 lines in total), `tests/`
(5 subtrees) and `docs/`.

The asymmetry starts here and never goes away. ECOGEM is a module of a model.
MARBL is a library with no model around it.

## 2. The trait space: the row's claim verified, with two corrections

`par_ecogem_plankton_file` is a namelist string, `ecogem_lib.f90:172-173`. The
file it names is read by `sub_init_populations`, `ecogem_data.f90:700-755`,
which counts the data lines, allocates `pft`, `diameter` and `random_n` to that
count, reads the rows, and sets `npmax=loc_n_elements` at `ecogem_data.f90:752`.
Every plankton-dimensioned array in the module is `ALLOCATABLE` and sized from
`npmax`, `ecogem_lib.f90:284-322`. **The population count is a runtime property
of a text file. The row's "file swap with no recompilation" is exact.**

The shipped ladder, counted between the `-START-OF-DATA-` and `-END-OF-DATA-`
tags of each file in `vendor/cgenie/genie-ecogem/data/input/`:

| file | populations | composition |
| --- | ---: | --- |
| `1P.eco` | 1 | one 19 um phytoplankton |
| `NPD.eco` | 1 | one 10 um phytoplankton |
| `NoDiat4ZP_PiEu.eco` | 8 | 2 Picoplankton, 2 Eukaryote, 4 Zooplankton |
| `9M.eco` | 9 | 9 Mixotroph, 0.2 to 1900 um |
| `3Diat4ZP_PiEu.eco` | 11 | the 8 above plus 3 Diatom |
| `8P8Z.eco` | 16 | 8 Phytoplankton and 8 Zooplankton, 0.6 to 1900 um |
| `8P7Z1F.eco` | 16 | 8 Phytoplankton, 7 Zooplankton, 1 Foram |
| `8P7Z3F.eco` | 18 | 8 Phytoplankton, 7 Zooplankton, 3 Foram |
| `32P32Z.eco` | 64 | 32 Phytoplankton and 32 Zooplankton |

**Correction 1: 0.6 to 1900 um is `8P8Z`'s spectrum, not the ladder's.** The row
attaches that range to the trait space in general. It belongs to one file.
`8P8Z.eco` is eight half-decade steps, 0.6, 1.9, 6, 19, 60, 190, 600, 1900 um,
mirrored exactly between the two functional types. `32P32Z.eco` is a different
and wider grid: phytoplankton from 0.10 to 1224.00 um, zooplankton from 0.15 to
2800.00 um, and the two are offset by one class rather than mirrored. Going from
sixteen populations to sixty-four therefore changes the size RANGE as well as
the resolution, which is a different experiment from refining one spectrum.

**Correction 2: the matched pair is matched in intent, not in count.**
`3Diat4ZP_PiEu.eco` has 11 rows against `NoDiat4ZP_PiEu.eco`'s 8; the three
extra rows are the diatoms, at 2, 20 and 200 um. Every other row is identical.
It is a diatoms-present against diatoms-absent pair, which is what the row calls
it, but the two configurations do not carry the same number of populations, so
the difference between them is a diatom guild ADDED rather than a diatom guild
switched off. That distinction is exactly the kind a pre-registered structural
sensitivity has to state.

The PFT name is matched case-insensitively, `ecogem_data.f90:329` calls
`lower_case`, then a chain of eleven string comparisons at
`ecogem_data.f90:330-419` sets six binary or near-binary traits per population:
`NO3up`, `Nfix`, `calcify`, `silicify`, `autotrophy`, `heterotrophy`, plus
`palatability` for some. The recognised names are prochlorococcus,
synechococcus, picoplankton, picoeukaryote, diatom, coccolithophore,
diazotroph, phytoplankton, eukaryote, zooplankton, mixotroph and foram. The
error message at `ecogem_data.f90:426` lists ten of those twelve, omitting
Eukaryote and Foram, both of which the code above it accepts and both of which
shipped `.eco` files use. The message is stale; the dispatch is not.

## 3. What the trait space actually DECLARES, and it is not the data file

This is the part OCN-4 needs and the row does not reach. The `.eco` file
declares only size and functional-type identity. Everything continuous is
generated from cell volume by allometry, in `sub_init_plankton`,
`ecogem_data.f90:311-697`:

    volume(:) = 1.0/6.0 * const_pi * diameter(:) ** 3     ecogem_data.f90:322
    qmin(iNitr,:) = qminN_a * volume(:) ** qminN_b        ecogem_data.f90:365
    ... 39 statements of the form coef_a * volume ** coef_b ...
    biosink(:)    = biosink_a * volume(:) ** biosink_b    ecogem_data.f90:539

`ini_ecogem_nml` declares 36 `_a` coefficients and their `_b` exponents,
`ecogem_lib.f90:62-152`. Maximum photosynthetic rate is the exception, a
three-parameter log-quadratic in volume at `ecogem_data.f90:466`, with five
per-PFT power-law overrides at `ecogem_data.f90:468-472`. The DOM partition
fractions are a two-asymptote sigmoid in DIAMETER rather than volume,
`ecogem_data.f90:565-566`.

**So the trait space Vesper would have to declare is a set of allometric
coefficient pairs, not a table of populations.** Adding a population costs one
line of text and no parameters. Changing the physiology costs a re-derivation of
36 coefficient pairs, every one of which is an Earth fit.

MARBL is the structural opposite. `autotroph_settings` is an array of a derived
type with **38 fields per autotroph**, `defaults/settings_latest.yaml:890-1272`:
sname, lname, temp_func_form_opt, Nfixer, imp_calcifier, exp_calcifier,
silicifier, is_carbon_limited, kFe, kCO2, kPO4, kDOP, kNO3, kNH4, kSiO3,
Qp_fixed, SiOpt, gQfe_max, gQfe_min, FeOpt, gQp_max, gQp_min, POpt, gQn_max,
gQn_min, NOpt, alphaPI_per_day, PCref_per_day, thetaN_max, loss_thres,
loss_thres2, temp_thres, mort_per_day, mort2_per_day, agg_rate_max,
agg_rate_min, loss_poc, Ea. Zooplankton carry 8 fields each
(`settings_latest.yaml:1274`); grazing relationships carry 13 each over a
`max_grazer_prey_cnt` by `zooplankton_cnt` array (`settings_latest.yaml:1349`).

`PFT_defaults` is a string with three valid values, `settings_latest.yaml:290-302`:
`CESM2`, `None`, `user-specified`. `autotroph_cnt` defaults to 1, is 3 under
`CESM2`, and is annotated `cannot change : PFT_defaults == 'CESM2'` and
`must set : PFT_defaults == 'user-specified'`, `settings_latest.yaml:852-863`.
The shipped `settings_latest+4p2z.yaml` proves the mechanism generalises:
`autotroph_cnt` 4, `zooplankton_cnt` 2 under a `4p2z` preset, lines 858-886 of
that file.

**So the row's "MARBL is locked at three" is too strong.** MARBL will take any
autotroph count. What it will not do is GENERATE the traits. Under
`user-specified` the default for every per-autotroph real is `1e34`
(`settings_latest.yaml:974` and following), a poison value, so an
N-autotroph configuration is a hand-written 38N + 8M + 13PM parameter table.
There is no size axis, no allometry, and no master trait: the identity of a PFT
is its parameter row.

**That is the real scientific difference, and it is sharper than the row's
version.** ECOGEM makes composition emerge from a low-dimensional declared trait
space; MARBL makes composition a consequence of a declared parameter table. On
the mechanism OCN-4 names as its reason for wanting trait-based ecology, ECOGEM
wins, and the reason is the allometric generator rather than the file swap.

## 4. State, restart and diagnostic cost

ECOGEM's prognostic state is `plankton(iomax+iChl, npmax, n_i, n_j, n_k)`,
allocated in `initialise_ecogem.f90`, where `iomax` is the number of elemental
quotas actually switched on. The quota switches are `nquota`, `pquota`,
`fquota`, `squota`, `chlquota` at `ecogem_lib.f90:42-47`, and the indices are
built by the `MERGE` chain at `initialise_ecogem.f90:100-106`, so unused
elements cost nothing.

The netCDF restart writes one 3-D field and one 2-D surface field per
(quota, population) pair, skipping chlorophyll for non-autotrophs,
`ecogem_data_netCDF.f90:101-118` and `131-145`. Nutrients are not in it; BIOGEM
owns those. Three points on the cost curve, at full C-N-P-Fe-Si plus
chlorophyll:

| configuration | quotas | 3-D restart fields |
| --- | ---: | ---: |
| `1P` | 6 | 6 |
| `8P8Z` | 6 | 5*16 + 8 = 88 |
| `32P32Z` | 6 | 5*64 + 32 = 352 |

The shipped `muffin.CBE.p_worjh2.BASESFeTDTL.FeMIP.SPIN` userconfig runs `8P8Z`
with `eg_nquota=.false.` and `eg_useNO3=.false.`, so C, P, Fe and chlorophyll
only: 3*16 + 8 = 56 3-D fields. Population count and quota count multiply, and
the quota count is the cheaper of the two to reduce.

MARBL's tracer set is enumerated in `defaults/settings_latest.yaml` under
`_tracer_list`: 17 unconditional base-biology tracers (PO4, NO3, SiO3, NH4, Fe,
Lig, O2, DIC, DIC_ALT_CO2, ALK, ALK_ALT_CO2, DOC, DON, DOP, DOPr, DONr, DOCr),
then per-autotroph Chl, C and Fe, plus Si for silicifiers, CaCO3 for calcifiers,
P when `lvariable_PtoC`, N when `lvariable_NtoC`; plus one C tracer per
zooplankton; plus optional abiotic-DIC and carbon-isotope modules.

**MARBL's restart is not MARBL's.** Tracers belong to the host. What MARBL
itself asks to have persisted is `marbl_saved_state`, and it is five fields at
most: PH_SURF, ABIO_PH_SURF, PH_SURF_ALT_CO2, PH_3D and PH_3D_ALT_CO2,
`src/marbl_saved_state_mod.F90:43-96`. That is an unusually clean restart
contract and it is the strongest single piece of evidence for the row's
engineering verdict.

## 5. The coupling shape, which is where the two are least comparable

ECOGEM's loop is `subroutine ecogem`, `ecogem.f90:3-32`. Seven arguments:
timestep length, the GENIE clock in milliseconds, incident shortwave, mixed
layer depth and sea-ice fraction on the (i,j) grid, the ocean tracer field, and
two output arrays. Its dimensions `n_i`, `n_j`, `n_k`, `n_ocn` and `n_sed` come
from `gem_cmn`, not from arguments. Initialisation takes GOLDSTEIN's vertical
grid and level-1 index directly, `initialise_ecogem.f90:5-30`, and
`check_egbg_compatible` at `ecogem_box.f90:517-568` stops the run if BIOGEM has
not selected an ocean tracer that ECOGEM's nutrient switches demand.

**ECOGEM cannot be lifted out of cGENIE.** Its host is not a parameter. What it
costs this project is whatever adopting cGENIE costs, and after that a namelist
block and a text file.

MARBL is the opposite shape and the burden is explicit. `marbl_domain_type`,
`src/marbl_interface_public_types.F90:54-66`, is `num_PAR_subcols`,
`num_elements_surface_flux`, `num_elements_interior_tendency`, `km`, `kmt`,
`zt`, `zw`, `delta_z`. There is no horizontal grid and no latitude anywhere in
`src/`. MARBL computes tendencies for columns and nothing else. The host owes it,
from `src/marbl_interface.F90:71-146`:

- the domain and the tracer index metadata;
- `tracers(:,:)` and `tracers_at_surface(:,:)`, which the host also time-steps;
- `bot_flux_to_tend(:)`;
- up to 16 named surface-flux forcings, `src/marbl_init_mod.F90:460-565`:
  u10_sqr, sss, sst, Ice Fraction, Dust Flux, Iron Flux, NOx Flux, NHy Flux,
  external C/P/Si Flux, Atmospheric Pressure, xco2, xco2_alt_co2, d13c, d14c;
- the interior-tendency forcings, `src/marbl_init_mod.F90:638-745`: Dust Flux,
  PAR Column Fraction, Surface Shortwave, Potential Temperature, Salinity,
  Pressure, Iron Sediment Flux, Iron Red Sediment Flux, Iron Vent Flux, O2
  Consumption Scale Factor, Particulate Remin Scale Factor, plus one restoring
  field per restored tracer;
- persistence of the five saved-state fields across restarts;
- **global reductions.** `glo_avg_fields_*` come out, `glo_avg_averages_*` go
  back in, and `set_global_scalars` at `src/marbl_interface.F90:1098-1114` feeds
  them to the bury-coefficient adjustment. A comment at
  `src/marbl_interface.F90:114` marks the running means as computed in the
  driver "for now". A column library that needs a global mean from its host is
  a real interface obligation, and it is the one most likely to be missed.

**Correction 3: the standalone driver is a test harness, not an offline model.**
The row and survey section 16b both read `tests/driver_src` as meaning MARBL
"runs without an ocean model against prescribed forcing", and infer that it can
be tried before an ocean host exists. Read at source, `tests/driver_src/marbl.F90`
dispatches on `testname` through ten cases, lines 304-573: init, init-twice,
gen_settings_file, request_diags, request_tracers, request_forcings,
request_restoring, bury_coeff, available_output, call_compute_subroutines,
get_put, marbl_utils. The one that computes anything,
`marbl_call_compute_subroutines_drv.F90`, reads a fixed NetCDF file of initial
conditions and forcing (`marbl_call_compute_subroutines_drv.F90:38`), calls
`surface_flux_compute` once at line 282 and `interior_tendency_compute` once per
column at line 342, and writes the results. **There is no time loop and no
tracer update.** It evaluates tendencies from a given state, once.

That does not destroy the sequencing argument, but it narrows it precisely. What
MARBL can be made to answer without a host is "how much does this tendency move
when I change this forcing", which is exactly the shape of a light-limitation
sensitivity test. What it cannot produce without a host is a spun-up state, an
equilibrium, or a productivity field. Any sentence that says MARBL "can be run"
should say which of those two it means.

## 6. The photon axis, re-checked, and it goes the other way

Survey section 12c concludes that MARBL is better on the K-star photon question
because "the conversion is the driver's responsibility and never a constant
inside the model", and the row records the axis as neutralised because `PARfrac`
is settable. The first half of that is wrong at source.

**Correction 4, and it inverts the verdict on the settability criterion.**

MARBL:

    real(r8), parameter :: f_qsw_par = 0.45_r8   ! PAR fraction
                                        src/marbl_settings_mod.F90:160
    PAR%interface(0,:) = f_qsw_par * ...
                                        src/marbl_interior_tendency_mod.F90:615

`f_qsw_par` is a Fortran `parameter`. It appears in the YAML settings files
nowhere, has no `put_setting` binding, and is used at exactly one site. The host
supplies total shortwave and an areal column fraction; **the shortwave-to-PAR
fraction itself is compile-time and is not reachable from any settings file.**

MARBL's water-column attenuation is likewise two hardcoded solar fits,
`src/marbl_interior_tendency_mod.F90:646-651`:

    if (WORK1(k) < 0.13224_r8) then
      PAR%KPARdz(k) = (0.0919_r8*unit_system%len2m)*(WORK1(k)**0.3536_r8)
    else
      PAR%KPARdz(k) = (0.1131_r8*unit_system%len2m)*(WORK1(k)**0.4562_r8)
    end if

Six literals, none in the YAML.

ECOGEM:

    PAR(:,:) = PARfrac * dum_egbg_fxsw(:,:)          ecogem.f90:129
    k_tot    = (k_w + k_chl*totchl)                  ecogem.f90:313, 320

`PARfrac`, `k_w` and `k_chl` are all in `ini_ecogem_nml`, `ecogem_lib.f90:106-110`,
with defaults 0.43, 0.04 and 0.03 in
`vendor/cgenie/genie-main/src/xml-config/xml/definition.xml:9031-9048`.

So on reachability the ranking is the reverse of the survey's: ECOGEM exposes
both the band fraction and both attenuation coefficients to a namelist; MARBL
hardcodes the band fraction and all six attenuation coefficients as Fortran
parameters.

MARBL is nonetheless better on STRUCTURE, and the two claims do not conflict.
Its attenuation is chlorophyll-dependent with a regime break, computed per layer
and per `delta_z`, against ECOGEM's single scalar pair. Its PAR is propagated
interface by interface with a floor, lines 660-678. But `KPARdz` has no
sub-column dimension, and `col_frac` is documented at
`src/marbl_interface_private_types.F90:24` as "column fraction occupied by each
sub-column". **`num_PAR_subcols` is an areal decomposition, not a spectral one.**
Both models are single-band below the surface, which is the defect survey
section 12c identifies in ECOGEM and it is present in MARBL too.

One constant nobody has recorded, in ECOGEM, on the same path:

    E0 = PARlocal/0.2174  ! convert from W m^-2 to muEin m^-2 s^-1
                                        ecogem_box.f90:298

1/0.2174 is 4.5998 umol photons per joule, the standard solar value. It is a
literal inside `photosynthesis`, not a namelist parameter. So ECOGEM carries TWO
solar photon constants on its light path and exposes only one of them. Survey
section 12a's measurement bounds what the buried one is worth: the mean photon
energy inside 400 to 700 nm moves by 2.3 percent between a 5772 K Planck and
k25v, against 15 percent for the band fraction. The survey's conclusion that
`PARfrac` being settable neutralises the axis therefore survives; what needed
correcting is that it is not the only constant there.

Two clear-water e-folding depths, computed from the constants above and stated
so the comparison is on the record: ECOGEM's `k_w` of 0.04 per metre gives
25.0 m; MARBL at its chlorophyll floor of 0.02 gives 43.4 m, falling to 22.3 m
at the regime break. Both are solar calibrations, and survey section 40e's
finding is that this star's flux below the surface layer is a smaller fraction
of the total. Neither model's clear-water constant can express the reason.

## 7. Fixed Earth relations that survive, and the stoichiometry claim does not hold

The row records stoichiometry as "a wash", on the grounds that MARBL exposes
`parm_Red_*` in YAML and BIOGEM exposes `par_bio_red_*` in its namelist, so
neither buries them. Both halves need correcting, in opposite directions.

**MARBL exposes three Redfield parameters, not eight, and derives seven.** The
module declares ten, `src/marbl_settings_mod.F90:300-311`. Three are settable:
`parm_Red_D_C_P` (112), `parm_Red_D_N_P` (16) and `parm_Red_Fe_C` (3.0e-6),
`src/marbl_settings_mod.F90:479-482` and `1291-1322`, and those three are the
only ones in `defaults/settings_latest.yaml`, lines 814, 823 and 839. The other
seven are computed at `src/marbl_settings_mod.F90:2158-2193` and printed as
derived, including `parm_Red_P_C_P = parm_Red_D_C_P` at line 2162, which asserts
that particulate and dissolved C:P are equal, and
`parm_Red_D_C_O2_diaz = parm_Red_D_C_P/150.0_r8` at line 2191, which carries a
bare literal. Those relations are fixed and none of them is settable. Three free parameters with a stated closure is a better
position than eight independent knobs, but it is a different claim.

**ECOGEM does bury Redfield, in the arithmetic that returns fluxes to BIOGEM.**
Four literals, none namelist-reachable:

    dum_egbg_sfcdiss(io_O2,:,:,:) = - 138.0 / 106.0 * nutrient_flux(iDIC,...)
                                        ecogem.f90:778
    dum_egbg_sfcdiss(io_ALK,:,:,:) = -16.0 * nutrient_flux(iPO4,...)
                                        ecogem.f90:783
    VCN(:) = up_inorg(iPO4,:) * 16.0
                                        ecogem_box.f90:311
    plankton(iChlo,jp,...) = chl2nmax / 6.625 * plankton(iCarb,jp,...)
                                        ecogem_data.f90:658

138/106 is the O2:C respiration ratio; the two 16.0 are N:P; 6.625 is 106/16.
The alkalinity line is the fallback branch taken whenever nitrate is not a
selected tracer, which is the configuration the shipped FeMIP userconfig runs.
A fifth, `up_inorg(iPO4,:) * 40.0` at `ecogem_box.f90:309`, sets the diazotroph
N:P, and its own code comment says it "probably should scale with nitrogen
fixation rate".

The irony is that ECOGEM's PLANKTON stoichiometry is emergent, through dynamic
`qmin`/`qmax` quotas per element per population, which is more flexible than
MARBL's; and MARBL's autotroph stoichiometry is also variable, with
`lvariable_PtoC` and `lvariable_NtoC` both defaulting true
(`defaults/settings_latest.yaml:359-373`) and `gQp_max`/`gQp_min`/`POpt` and
their N and Fe counterparts in the autotroph table. **The difference is not in
the ecology at all. It is that ECOGEM's COUPLING arithmetic reverts to fixed
Earth ratios wherever the quota it needs is switched off, and there is no
warning when it does.**

Other fixed relations found and worth having on the record:

- ECOGEM's grid area uses `const_rEarth` directly, `ecogem_data.f90:884`, which
  reads `2.0*const_pi*(const_rEarth**2)*(1.0/n_i)*(sv(j)-sv(j-1))`.
  `const_rEarth` is `6.37E+06`, `vendor/cgenie/genie-main/src/fortran/cmngem/gem_cmn.f90:802`.
  It is used only for the time-series site diagnostic, not for the tendencies.
- ECOGEM's ballast partition, `ecogem.f90:817-822`, uses literal 0.085, 0.025
  and 0.0 where the commented-out line directly above it used the namelist
  parameters `par_bio_remin_kc`, `par_bio_remin_ko` and `par_bio_remin_kl`. The
  namelist route was replaced by literals and the parameters left in place.
- MARBL's `dust%diss = 400.0_r8` and `P_iron%diss = 600.0_r8`,
  `src/marbl_interior_tendency_mod.F90:974` and `980`, are dissolution length
  scales in metres set in Fortran, with no YAML entry, unlike `parm_POC_diss`,
  `parm_SiO2_diss` and `parm_CaCO3_diss` beside them. **MARBL's settings schema
  is not complete**, which is worth knowing before anyone treats it as the
  authoritative parameter inventory.
- MARBL's temperature functions offer `q_10`, `arrhenius` and `power` per PFT,
  `defaults/settings_latest.yaml:917-925`. ECOGEM's `ctrl_tdep_form` offers an
  exponential form, Eppley and MEDUSA, `ecogem_box.f90:215-232`, with
  `temp_T0` a reference temperature defaulting to 20.0 degrees Celsius,
  `definition.xml:9197`. Both are Earth calibrations and both are switchable.

## 8. Absolute time, and it is compile-time in both

OCN-4 asks for the absolute-time dependencies. Both models state rates per day
and convert with a hardcoded seconds-per-day.

ECOGEM: `real,parameter :: pday = 86400.0`, declared **three separate times**, at
`ecogem_data.f90:320`, `ecogem.f90:49` and `ecogem.f90:886`. Every rate
parameter in `ini_ecogem_nml` is divided by it at `ecogem_data.f90:644-650`, and
diagnostics are multiplied back by it at `ecogem.f90:729` and `907-916`. There is
no year in ECOGEM's own source; it takes `conv_yr_s` from `gem_cmn`
(`gem_cmn.f90:587`) and there is one further literal, the 48 in
`loc_yr = real(dum_genie_clock)/(48.0*1000.0*conv_yr_s)`, `ecogem.f90:985`, which
is BIOGEM's timesteps per year restated inside ECOGEM. That is a further instance
of the pattern survey section 46b records for cGENIE as a whole.

MARBL: `spd = 86400.0` and `dpy = 365.0` with `spy = dpy*spd`,
`src/marbl_constants_mod.F90:38-47`, in one module, used by the `_per_day` and
`_per_year` named parameters and by `c14_lambda` and the DOM remineralisation
rates (`src/marbl_settings_mod.F90:164-166`, "1/15yr"). One statement of the day,
one of the year, both compile-time.

Neither is namelist-reachable. MARBL's is centralised and named; ECOGEM's day is
triplicated and its year is somebody else's. On this axis MARBL is cleaner and
the difference is one of maintenance, not of capability, since both need a source
edit.

## 9. Sinking and gravity, where the two are structurally different and ECOGEM's is dead

Neither model contains a gravitational constant. `grep -rn 'grav'` over
`references/marbl/src/*.F90` returns nothing, and ECOGEM has none either. Both
express particle export empirically, but in DIFFERENT DIMENSIONS, and that
matters for a planet with a different surface gravity.

MARBL uses e-folding LENGTHS. `decay_POC_E = exp(-dz_loc / poc_diss)`,
`src/marbl_interior_tendency_mod.F90:3085`, applied at line 3190, with
`POC%diss = parm_POC_diss` at line 944, default 100 m
(`src/marbl_settings_mod.F90:451`, expressed as `100.0e2` cm and converted at
line 485). CaCO3 500 m, SiO2 650 m. A length scale is a sinking speed divided by
a remineralisation rate, and MARBL exposes neither factor separately, so there is
no defensible way to rescale it for a different gravity from inside the model.

ECOGEM declares a sinking VELOCITY, `biosink(:) = biosink_a * volume(:) **
biosink_b` at `ecogem_data.f90:539`, converted to per-second at line 649. A
velocity does have a first-principles gravity scaling in the Stokes regime.

**But ECOGEM's is dead code.** A whole-tree grep for `biosink` returns eleven
hits and every one is a declaration, an allocation, the namelist, the
computation, the unit conversion, or a write to a diagnostic file. It is never
read by any tendency. `ecogem.f90:363` carries the commented-out call site,
`!call sinking`. The shipped defaults are `biosink_a = 0.0` and `biosink_b = 0.0`,
`definition.xml:9135-9140`, so the parameter would be zero even if the path
existed.

The finding this establishes is cleaner than either half: **in the ECOGEM tier,
sinking is not ECOGEM's.** ECOGEM partitions mortality and messy feeding into
dissolved and particulate pools and hands the particulate flux to BIOGEM
(`ecogem.f90:760-772`); BIOGEM's remineralisation profile does the vertical
transfer. Any gravity question about export in this tier is a BIOGEM question,
and `biosink` is a parameter that looks like an answer and is not one.

## 10. Three defects in the vendored ECOGEM, found while reading

Recorded because they change what a reader can trust in the tree, not because
anything here proposes to fix them.

**The grazing DOM partition uses the mortality parameter.**

    beta_graz(:) = beta_graz_a - (beta_graz_a-beta_graz_b) / (1.0+beta_mort_c/diameter(:))
                                        ecogem_data.f90:565
    beta_mort(:) = beta_mort_a - (beta_mort_a-beta_mort_b) / (1.0+beta_mort_c/diameter(:))
                                        ecogem_data.f90:566

`beta_graz_c` is declared at `ecogem_lib.f90:143`, exposed in the namelist at
`ecogem_lib.f90:148`, printed to the log at `ecogem_data.f90:135`, and never read
by anything. The grazing line uses `beta_mort_c`. Every shipped userconfig sets
both to 100.0, so the defect is latent rather than live in the configurations
that exist. A configuration that separated them would silently get the mortality
value in both places.

The same two lines carry a units question that the source cannot settle. Every
userconfig comment reads "Size at 50:50 partition (default = 100 um^3)", a
volume, and the code divides by `diameter(:)`, which is in um. Either the comment
is wrong or the intended variable was `volume`. The two readings differ by the
cube of a length, so this is not a cosmetic discrepancy. **Bracketed rather than
guessed**: the source is consistent with either, and resolving it needs Ward et
al. (2018) section-level reading, which was not done here.

**The shutdown output is jumped over.** `end_ecogem.f90:18` is an unconditional
`goto 101`, and `101 continue` is at line 35. Everything between, the 2-D netCDF
write and the disabled 3-D block, is unreachable. The subroutine's own header
comment at line 3 reads "END GEMlite".

**The shipped default plankton file does not exist.**
`definition.xml:9234-9235` sets `par_ecogem_plankton_file` to `one_plankton.eco`.
A `find` over `vendor/cgenie` returns eleven `.eco` files and none of them is
that one. `sub_init_populations` stops the run when the file yields zero rows,
`ecogem_data.f90:714-720`. Every working userconfig overrides the key. Since the
file swap IS the trait-space interface, a default that names a missing file is
worth knowing about before anyone treats the shipped default as a starting
configuration.

## 11. What this establishes, and where the row's conclusion lands

The row's central claim survives: **the comparison splits scientific against
engineering, and the evidence does not support forcing it to a winner.** But
both halves are different from the row's version of them.

On the science, ECOGEM's advantage is larger than the row states and rests on a
different mechanism. It is not the file swap, which is a convenience; it is the
allometric generator, which means a Vesper trait space is declarable as roughly
36 coefficient pairs plus a size grid, against MARBL's 38-parameter table per
population with no generator and a poison default. Sixty-four populations is the
top of the ladder, and reaching it also widens the size range.

On the engineering, MARBL's advantage is real, and the restart contract is its
strongest single piece of evidence: five persisted fields, with tracers owned by
the host. But its driver-agnosticism is the agnosticism of a LIBRARY. Its host
obligations are large and include a global reduction, and its standalone
executable evaluates tendencies once rather than integrating anything.

And the two are not independent choices. ECOGEM cannot be extracted from cGENIE:
its dimensions come from `gem_cmn`, its grid from GOLDSTEIN and its nutrient
supply from BIOGEM, which also owns the sinking that ECOGEM's own sinking
parameter pretends to. **So the ecosystem tier is entangled with the circulation
host in a way the row's framing does not show.** If cGENIE is adopted, ECOGEM
costs a namelist block and a text file. If it is not, ECOGEM is unavailable at
any price, and MARBL's cost is not the interface but the tracer transport model
underneath it, which this project does not have.

Three of the row's supporting claims should not be carried forward as stated:
MARBL is not locked at three autotrophs; MARBL's PAR fraction is not the
driver's responsibility but a Fortran parameter, which inverts the photon axis
on reachability while leaving the survey's verdict on magnitude intact; and
stoichiometry is not a wash, because ECOGEM buries four Redfield literals in its
coupling arithmetic and MARBL derives seven of its ten ratios from three.

## 12. What this reading did NOT establish

- **Neither model was compiled or run.** Every claim above is from source and
  from the shipped configuration files. Nothing here has been observed to
  execute.
- **The Si cycle is section 13**, added after this list was written. It reads
  both tiers' silicon paths at source and reads Naidoo-Bagwell et al. (2024),
  which was `held` when this section was first written.
- **The magnitude of the trait-space re-derivation is not bracketed.** That 36
  coefficient pairs are Earth fits is established; what any of them would be on
  Vesper is not, and nothing here estimates it.
- **MARBL's diagnostics inventory was not read.** `defaults/diagnostics_latest.yaml`
  is 2843 lines and only its existence is used above.
- **BIOGEM's remineralisation profile was read at line level only on the OPAL
  path**, in section 13c. Section 9's redirection of the gravity question to
  BIOGEM is answered there for opal and remains a redirection for the organic
  carbon and carbonate paths.
- **The `beta_graz_c` units question is open**, as stated, and is bracketed
  rather than resolved.

## 13. The silicon axis, closed

Section 12 recorded this as the one axis the reading had not covered, and the
row asks for the required C-N-P-S-Fe-Si inputs of each tier. Read at source
2026-08-24, with `references/naidoobagwell2024-ecogenie-diatom-extension.pdf`
read for the parts the source cannot settle.

### 13a. ECOGEM: silicon is namelist-reachable and gated three ways

Silicon is compiled in unconditionally; there is no preprocessor guard anywhere
in `vendor/cgenie/genie-ecogem/src/fortran/`. It is gated at runtime by two
namelist logicals and one per-population real: `squota`, the dynamic silicon
quota (`ecogem_lib.f90:45`, namelist at `:47`, default false); `useSiO2`,
silicate as a resource (`:38-39`, default false); and `silicify(npmax)`,
allocated at `initialise_ecogem.f90:229`.

**The silicifier is selected by a word.** `ecogem_data.f90:364` sets
`silicify(jp) = 1.0` only where the functional type name is `diatom`, matched
case-insensitively; the other eleven recognised names set zero.

**The `.eco` files carry no silicon column at all.** Every shipped file has three
columns: functional type, cell diameter in micrometres, replicate count. The
silicon traits come entirely from the name in column one and the diameter in
column two, through five allometric `_a`/`_b` pairs declared at
`ecogem_lib.f90:94-101` and applied at `ecogem_data.f90:522-528`, all of them in
`ini_ecogem_nml`.

**One of the nine shipped plankton files carries a silicifier row**,
`3Diat4ZP_PiEu.eco:6-8`, at 2, 20 and 200 micrometres, and one base config in the
whole tree combines ECOGEM with the silicon tracers. Section 2's correction to
the row stands and sharpens: the diatom and no-diatom pair is a guild ADDED, and
the shipped pair of user configs differs by exactly one line, the plankton file
name, so in the no-diatom arm silicon is a fully enabled and completely inert
tracer.

Four things the source settles that a structural trace did not:

- **There is no silicon half-saturation parameter; it is derived.** The uptake
  law at `ecogem_box.f90:165-166` has Michaelis-Menten equivalent
  `k = vmax/affinity`, and the model prints exactly that at
  `ecogem_data.f90:112`.
- **The Si:C quota is size-INDEPENDENT while the kinetics are not.**
  `qminSi_b` and `qmaxSi_b` default to zero and no user config in the tree sets
  either, while the uptake rate and affinity exponents are both nonzero. A 2 and
  a 200 micrometre diatom have the same Si:C and different kinetics.
- **Silicon limitation is linear, not Droop** (`ecogem_box.f90:113`), and enters
  the von Liebig minimum only for silicifiers (`:122`). Zooplankton assimilate no
  silicon (`ecogem.f90:402-404`).
- **`kexcSi` is dead code**, computed at `ecogem_data.f90:528` and read by
  nothing, the same class as `biosink` in section 9. Its default is zero so this
  is latent, but a Vesper retune that set a silicon excretion rate would silently
  get nothing. `par_diatom_vmax_mod` is a second dead knob, namelisted at
  `ecogem_lib.f90:236-237` and never read, and it is the obvious one a diatom
  retune would reach for.

### 13b. ECOGEM: one conservation requirement is enforced and one is not

`ecogem_box.f90:553-566` hard-stops if `ocn_select(io_SiO2)` is off. **Nothing
checks `sed_select(is_opal)`.** ECOGEM writes the opal flux unconditionally on
`squota` at `ecogem.f90:770`, and BIOGEM copies only SELECTED sediment tracers
(`biogem.f90:2388-2402`). With opal unselected, silicon uptake is a permanent
removal with no return path and nothing in the tree catches it.

**Dissolved organic silicon does not exist.** There is no `io_DOM_Si` tracer.
ECOGEM computes the size-dependent dissolved and particulate partition for the
silicon index anyway, its own comment at `ecogem.f90:542` saying so, and then
recombines both halves into opal at `:770`. So section 10's open units question
about the partition coefficients has no effect on the silicon path at all.

Two latent numerical defects, neither checked by `check_egbg_compatible`:
`squota` true with `useSiO2` false indexes `vmax(0,:)` at
`ecogem_data.f90:526`, and `useSiO2` true with `squota` false indexes `qreg(0,:)`
at `ecogem_box.f90:180`. Both are out of bounds. And the "Population inviable!"
guard at `ecogem_data.f90:525` is
`maxval(qmin(iSili,:)/qmax(iSili,:)) > 1.0`, which is zero over zero for every
non-silicifier and therefore for eight of the eleven populations in the shipped
configuration. Under the optimised build that is a NaN and the comparison is
false, so **the guard has no failing state**, which is the class this project's
own conventions forbid.

### 13c. Where the silicon gravity question routes, and it answers section 9

Section 9 redirected the gravity question from ECOGEM to BIOGEM and called that a
redirection rather than an answer. **On the opal path BIOGEM answers it.** Under
the kinetic dissolution branch the layer residence time is
`dz / sinkingrate_reaction` (`biogem_box.f90:3425`) and dissolution is
`1 - exp(-dt_reaction * par_bio_remin_opal_K * f(T,u))` (`:3482-3489`). **The
sinking SPEED and the dissolution RATE are two separate namelist-reachable
numbers and only the speed is a gravity question**, which is the seam the
survey's rate-time-base finding says must exist.

MARBL's 650 m has no such seam. It is the ratio of the two, collapsed, and there
is no defensible way to rescale it from inside the model.

The live speed in the published configuration is
`par_bio_remin_sinkingrate_reaction` at 49.689 m/d, not the 125.0 the
configuration also carries: `par_bio_remin_sinkingrate` is read only by a
backwards-compatibility branch (`biogem_data.f90:608-612`) that the config
disables. **Bracketed** on the same terms as
`notes/audits/ocean-tier-implicit-earth.md` finding 14: a Stokes scaling gives
about 64.9 m/d at 1.306 gravity, the aggregate regime gives less, the model
carries one scalar speed for all particles under its default remineralisation
function, and the settling regime is not determined here.

### 13d. MARBL: silicon is one tracer, and opal is not a tracer at all

`SiO3` is unconditional. A per-autotroph silicon tracer is created only where
`silicifier` is true
(`references/marbl/src/marbl_interface_private_types.F90:1551-1554`), and every
downstream test is on that index, so turning the flag off removes the whole opal
pathway with no other signal. **Opal itself is a
`column_sinking_particle_type`** (`:244`, constructed at `:919`), diagnosed
within a single column call, so it has no prognostic inventory and is implicitly
assumed to reach the seafloor inside one host timestep.

Three per-autotroph silicon fields exist and all are settings-reachable, with the
uptake half-saturation and the optimal silicate carrying the `1e34` poison
default under `user-specified` so that leaving them unset aborts. **`silicifier`
does not**: it defaults to false
(`references/marbl/src/marbl_pft_mod.F90:247`), so a user-specified plankton set
that forgets it gets zero silicifiers, no opal cycle, and no poison value to trip
on. That is the one place MARBL's engineering advantage does not hold.

**There is no independent silicon photosynthesis rate.** Opal production is the
carbon fixation rate times the silicon quota
(`marbl_interior_tendency_mod.F90:1824-1826`), so it inherits every constant on
the carbon fixation path and the diatom's light curve propagates into the opal
flux one for one. The transfer to opal uses a grazing remineralisation fraction,
`f_graze_si_remin = 0.22`, a bare Fortran `parameter` at
`marbl_settings_mod.F90:141` with no registration and no YAML entry, and the
comment three lines above the code that uses it says sixty per cent.

Burial is a step function at `marbl_interior_tendency_mod.F90:3492-3505`, 0.2 of
the flux above a threshold and 0.04 below, attributed in a comment to Ragueneau
et al. (2000); all three numbers are unnamed inline literals. The hard ballast
subclass decays on an unnamed `4.0e6` cm at `:3044`. Dust-derived silicate is
computed inside MARBL from two inline literals at
`marbl_surface_flux_mod.F90:434-435`, neither of them tunable.

**MARBL closes the column silicon budget and not the global one.**
`marbl_diagnostics_mod.F90:4480-4533` sums every silicon tendency plus the
sediment loss and aborts above a threshold, which is a real conservation check
and it is loud. The global closure is the bury-coefficient Newton step, every
input to which is a global average MARBL cannot compute. **Without a host
performing that reduction, MARBL's silicon cycle has no closure at all.**

And the standalone driver cannot supply one. Confirmed at source: twelve test
cases dispatched from one `select case`, no timestep and no step count in its
namelist, the compute routines called once per column, and the tracer array read
once and never written back. For silicon that means it can emit production,
remineralisation, formation and burial as instantaneous rates at one prescribed
state, and cannot produce an inventory, an annual opal production, a
steady-state dissolution profile, or any evidence that the cycle converges.

### 13e. What the paper settles that the source cannot

`references/INDEX.md` recorded Naidoo-Bagwell et al. (2024) as `held`. It is now
read, and five things in it are not readable from the tree.

- **The published EcoGEnIE 1.1 silicon cycle has no source and no sink, by
  declaration.** Section 2.2 states that the paper does not calculate opal
  preservation and imposes a benthic closure instead, and the base config has
  both the sediment and weathering components off. **The whole inventory is one
  configuration line.** Adopting that configuration adopts a closed silicon
  inventory.
- **The dissolution law's provenance.** The seven literals at
  `biogem_box.f90:3482-3489` carry no citation in the source; the paper names
  Ridgwell et al. (2002) and a sediment-trap evaluation.
- **The tuned parameters' ranges and their selection**, in its Table 2, which
  exist nowhere in the tree: a 550-member Latin hypercube scored against Earth
  ocean fields, and the accepted values reconcile exactly with the shipped user
  config.
- **The acceptance gate is an Earth inventory**: a global opal export compared
  against a published Earth range, used as a selection criterion.
- **Which parameter controls the silicifier guild's size structure**: the
  silicon uptake rate exponent, its section 5.2, which no reading of the source
  would surface.

**And one thing it does not settle, which matters here.** Its Table 2 attaches
two of the four live silicon parameter families to references that do not appear
in its own reference list, checked two ways. **So the literature basis for half
of ECOGEM's live silicon parameters cannot be followed from the paper**, and by
this project's own rule an undocumented value is a test case rather than an
anchor. Those tested ranges are not calibration targets.

### 13f. The comparison, and what would have to be true

| | ECOGEM | MARBL |
| --- | --- | --- |
| dissolved Si tracer | BIOGEM's, enforced | its own, always present |
| particulate Si | BIOGEM's opal and its refractory fraction; **not enforced** | none; diagnostic within a column |
| silicifier identity | a word in a text file | a per-autotroph logical defaulting to false |
| Si trait parameters | four live allometric pairs, all namelist | two per-autotroph plus three global, all settings-reachable |
| vertical transfer | a sinking SPEED times a dissolution RATE, separately reachable | one fixed e-folding LENGTH |
| Si source | weathering component, or nothing | host-supplied flux, or nothing |
| Si sink | sediment diagenesis, or nothing | its own burial step function plus a host global reduction |
| global reduction required | no | **yes** |
| what it answers without a host | nothing; it cannot leave cGENIE | one instantaneous tendency |

**On silicon the reachability comparison goes the same way as the photon axis,
and for the same reason.** Every live silicon parameter in ECOGEM is
namelist-reachable and its buried Earth content is one layer down in BIOGEM;
MARBL's grazing remineralisation fraction, burial step function, hard ballast
length and dust silicate pair have no settings entry at all. The one structural
advantage MARBL keeps is that its silicon budget check is real and loud, where
ECOGEM's missing sediment-tracer guard leaks silently.

Five preconditions before the silicon axis could be called settled:

1. **A silicon inventory would have to be DECLARED for Vesper.** In the ECOGEM
   tier it is literally one number and with weathering and burial off it is
   conserved exactly. Nothing in either model derives it, and it is the largest
   single piece of Earth on the axis.
2. **The half-saturations would have to be re-derived against that inventory
   rather than carried.** MARBL's optimal-silicate threshold sits at the edge of
   Earth's surface regime: raise the inventory and the entire variable Si:C
   mechanism switches off silently, lower it and the quota clamps to its minimum
   globally. The same argument applies to the uptake half-saturation, which
   decides whether diatoms are silicon-limited at all.
3. **The sinking speed would have to be settled**, and only cGENIE's tier has
   the seam to apply it at. Section 13c.
4. **The diatom light curve would have to move with the spectrum**, because opal
   production is the carbon fixation rate times a stoichiometric ratio in both
   models and neither has a silicon-independent light path. Section 6's factor
   applies unchanged, and ECOGEM makes it worse by confining all uptake to a
   single surface box chosen for a solar photic zone.
5. **Whether Vesper's ocean has a silicon source and sink at all is a design
   decision, not a model setting**, and the published configuration answers it
   with "neither".
