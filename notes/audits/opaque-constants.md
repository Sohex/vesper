# Opaque constants in the ExoPlaSim fork

Audited 2026-08-24. The question: where does a number reach the climate model
without the code, its comments, or any document in this tree being able to say
why it is that number? `inherited-earth-constants.md` asked a narrower version of
it over the radiation's surface and cloud terms, the orbit block and the flux
derivation. This one covers the whole fork -- all 46 files changed since the
subtree was added, plus the driver component in `exoplasim/` -- and separates
three things that had been running together: a constant nobody chose, a constant
whose derivation exists somewhere a reader of the code cannot reach, and a
constant whose documentation no longer describes it.

**The fork's own additions are the cleanest code in the tree.** The CO2 line-list
fit, the CH4/N2O band triples, `zovl`, `zvol2cm`, `zepsc`, the SHTns algorithm
selectors and the Gaussian-weight rewrite all carry their derivation in the code
plus a note or a JSON artifact; `legmod.f90` has no unexplained literal left in
it. Every finding below is therefore in inherited physics, in the instruments
that judge the model, or in values hand-transcribed from a generator.

**The two largest findings are not opaque constants at all but defects the audit
walked into**: the derived hyperdiffusion reaches one model level of ten, and
every damping timescale is applied over 1.25 times its intended interval. Both
are recorded here because both were found by asking where a number comes from.

---

## 1. The derived hyperdiffusion reaches one model level of ten

`run_exoplasim.py:1031-1034` and `stability_probe.py:177-181` write `TDISSD`,
`TDISSZ`, `TDISST`, `TDISSQ` and `NDEL` as bare scalars. Fortran namelist
assignment of a scalar to an array sets element 1 only. The model's own echo,
from `run_2b20e3324bb0/MOST_DIAG.00002`:

    NDEL=4          , 9*2          ,
    TDISSD= 0.4457     , 9*0.2000     ,
    TDISSZ= 2.2285     , 9*1.1000     ,
    TDISST= 5.6449     , 9*5.6000     ,
    TDISSQ= 0.7421     , 9*0.1000     ,

Levels 2 to 10 run the compiled T21 module defaults from `plasimmod.f90:877-880`.
That is the state `world-ys9` exists to end, still live after the derivation was
written.

**At T42 the array is mixed-unit.** There the `if(NTRU==42)` block at
`plasim.f90:1443-1450` has already filled levels 2 to 10 in SECONDS, and the
driver then overwrites level 1 in days. From `run_953ee807d32f`:

    TDISST=  2.8224     , 9*65664.0
    TDISSD=  0.2229     , 9*5184.0

`dayseccheck` (`plasim.f90:1643-1656`) discriminates on `maxval`, sees 65664
against a `deltsec` of 2700, and converts nothing. No "assuming [days]" line
appears in that run's diag. Level 1's temperature timescale is therefore read as
2.82 seconds against level 2 to 10's 65664: nondimensionally `tdisst(1)` = 4871
against 0.209, of order 23,000 times the damping on the top model level.

ExoPlaSim's own wrapper does this correctly --
`__init__.py:2804` writes `"%d*%f"%(self.layers,qdiffusion)` -- and the project's
scripts bypass it to edit the namelist key directly.

This also puts a second explanation under a result already on record.
`resolution-tuned-parameters.md` reports that applying the derived damping at T42
for two orbits did not move the KE spectrum, against a change in `tau_vorticity`
of 3.7x. The analytic filter-versus-diffusion comparison stands on its own; the
experiment that corroborated it did not apply what it recorded.

## 2. Every damping timescale is applied 1/rotspd too long

`plasim.f90:1783, 1789, 1799, 1804, 1809, 1814` and `1628` nondimensionalise with
`day_24hr` = 86400.0 (`plasimmod.f90:961`), while the nondimensional time unit is
`1/ww = sidereal_day/TWOPI` (`plasim.f90:1514`). The per-step damping works out to
`2*day_24hr*sak/(ntspd*tau)` against a correct `2*deltsec*sak/tau`, so the ratio
is `mtspd/ntspd` = 32/40 -- **exactly `rotspd`**, hence exactly 1.0 on Earth and
0.8 here. The applied timescale is 1.25 times the namelist's.

Verified against the run's own diag: `new maxval(tdisst) = 487719.36` = 5.6449 x
86400, and 86400/(TWOPI x 487719.36) = 0.0282 where 108000/(TWOPI x 487719.36) =
0.0353.

Applies identically to `restim`, `tfrc` and `dampsp`. Bounded on the diffusion by
`world-6on`'s finding that the filter dominates it by 450 to 1600 times; not
bounded on `tfrc`.

The upstream comment at `plasim.f90:1510-1512` calling `day_24hr` "purely used for
time unit conversions", and the unresolved `! day_24hr = 86400.0 ! WHY DOESN'T
THIS WORK???` at `:1502`, are therefore a false premise rather than a curiosity.
`day_24hr` cancels out of the timestep, but it is the denominator of this
nondimensionalisation and the unit that `tfrc`, `tdiss*`, `restim` and `dampsp`
are entered in. Changing it changes what the model integrates.

## 3. `t0 = 250.0 K` on every level, and the energy sink scales with it

`plasimmod.f90:882` declares `t0(NLEV) = 250.0`, the semi-implicit reference
temperature, isothermal on all ten levels. It is a `plasim_nl` key and no run
sets it: every diag echoes `T0= 10*250.0`. Nothing in `notes/`, `docs/src/` or
`bd` asks why 250, though `resolution-tuned-parameters.md` cites it a dozen times
as a given.

It ranks where it does because the project's own measurement makes the core
energy sink proportional to it. `Phi0_k = Sum_j g(j,k) t0` is built from `t0`
alone, and that note's finding is that `Phi0` weights the difference by up to 3.8
times `t0`, giving `Cimp - Ct = -0.9588 W/m2` against a sink of `-0.7938`. **The
defect the fork's energy fixer exists to mask scales with a number nobody
chose**, and 250 K is Earth's mass-weighted tropospheric mean rather than this
atmosphere's. It also sets the explicit gravity-wave limit `nconvtime` computes
(`plasim.f90:700`) and the blow-up floor at `plasim.f90:1338`.

## 4. The transform-equivalence gate certifies a filter the model stopped using

`verify_shtns_equivalence.f90:72` sets `nfilterexp = 8`, under the comment "as the
beds run it" and inside a block header asserting "THE CONFIGURATION UNDER TEST IS
THE ONE THE MODEL RUNS". `config/planet.yaml:388` declares `filter_power: 16`,
derived -- gamma >= 14.78 is required at T21 -- with `check_consistency.py`
recomputing it and failing if it is not met. It moved from 8 on 2026-08-23.
`filter_kappa` at line 71 still matches; only gamma is stale.

The direction matters. `exp(-8x^8)` damps the mid-band far harder than
`exp(-8x^16)`: 0.969 against 0.9999 at `n/NTRU` = 0.5, and the mid-to-high band is
where `shtns-viability.md`'s analysis arms carry their largest residual. **The
gate measures a smaller transform error than the shipped configuration would
produce.**

## 5. The CO2 shortwave fit exists in three places with two values

`radmod.f90:2262-2265` runs the HITRAN2020 refit under PHYS-10:

    zca1=3.0658E-4  zcb1=20.376  zca2=3.8193E-3  zcb2=3.9672E-3

`patches/exoplasim-3.4.2-co2-shortwave.patch:193-196` and
`analysis/shortwave_band_weights.json` (`co2.closed_form_fit`) both still carry
the Howard fit, `3.8265E-4, 44.539, 2.2325E-3, 5.8954E-3`, which is +7.7% in
absorptance at this planet's own CO2 path. `notes/shortwave-co2.md:163-165`
asserts the two agree, and `shortwave_band_weights.py:986-988` writes that claim
into every regenerated JSON, where it is true of the patch and false of the model.

Two consequences, in opposite directions: re-running the generator silently
reasserts the superseded fit as the stored artifact, and rebuilding from the patch
stack restores a CO2 shortwave term 7.2% too strong. `patches/README.md:20-22`
requires a superseded patch to say so; it does not, and it is upstream PR #62.

## 6. Surface longwave emissivity

`radmod.f90:3237`:

    zeps(:)=dls(:)+0.98*(1.-dls(:))

Land emissivity is exactly 1.0 by construction; ocean and sea ice are 0.98. No
comment, no unit, no citation. A grep for "emissivity" across `notes/`,
`exoplasim/notes/`, `docs/src/` and `bd search` returns nothing: **this is the
only large radiative literal in the module with zero coverage anywhere in the
tree.**

It is used consistently -- `:3238` for emission and `:3547` for the `(1-eps)`
reflection of downward longwave -- so it is not an energy-conservation defect. It
is an unexamined number. Its cost is `(1-eps)*(sigma*Ts^4 - LW_down)`, which at
`eps` = 1.0 over land is identically zero. Measured broadband emissivity of bare
rock and desert sand is 0.90 to 0.95, so on a world whose land-albedo work turns
on a large barren fraction that is roughly 3.5 to 7 W/m2 of overstated land
emission over the driest surfaces, one-signed. Ocean at 0.98 against seawater's
0.985 to 0.99 is roughly 0.5 to 1.5 W/m2.

## 7. The land roughness field is anchored through Earth's gravity

`build_surface_roughness.py:19, 221, 244` hardcodes 141.6 m as "this planet's
lowest level", three times, with no config key and no derivation.
`hydrography/scripts/carve_verdict.py:217` derives the same quantity as
`(GASCON * t_air / gravity) * ln(1/SIGMA_LOWEST)`.

Measured 2026-08-24, with `GASCON` 287.017 and `SIGMA_LOWEST` 0.9828:

| gravity | 279 K | 288 K |
| --- | --- | --- |
| 12.81 (this world) | 108.46 m | 111.95 m |
| 9.80665 (Earth) | 141.67 m | 146.24 m |

The hardcoded value is the Earth one, high by exactly the gravity ratio 1.306.

The builder bisects the free orographic coefficient so that the area-weighted land
mean of the derived `z0` field hits a target, through `ce = 0.16/ln(141.6/z0)^2`,
so the anchor is 14% low: `ce` at the 2.0 m default is 0.00882 as coded against
0.01004 correct. The "7.7x too much exchange" claim in the script's own docstring
and in `lake-representation.md:208-214` is really 8.4x. This reaches land
sensible heat and land evaporation, and therefore the carve criterion's numerator.
Separately, the `0.16` is an unstated von Karman squared.

The anchor itself is the second half: the target is ExoPlaSim's own `dz0land` =
2.0 m, whose only documentation in the model source is the comment "roughness
length land". `config-rationale.md:764-766` records the decision to anchor there
so the field "redistributes roughness without moving the global value the model
was tuned against", and neither that nor `steps.md:40-41` establishes that 2.0 m
is a tuned value rather than an unexamined default.

## 8. Cloud liquid water, and therefore cloud optical depth

`rainmod.f90:2005-2009`:

    zzh(:)=700.*ALOG(1.+dqvi(:))
    dql(:,jlev)=0.00021*EXP(-zzf(:,jlev)/zzh(:))*gascon*dt(:,jlev)/(sigma(jlev)*dp(:))

Two literals, no unit, no citation beyond a comment pointing at "CCM3
description" fifteen lines earlier. `dql` is the only thing that sets cloud
optical depth: `radmod.f90:2465-2466` builds `zlwp = 1000*dql*dp/ga*dsigma` and
then `ztau = 2.0*ALOG10(zlwp+1.5)**3.9`, and `radmod.f90:3270` uses it again for
the longwave.

`700.` is a length in Earth metres, the e-folding depth of cloud water, and it is
the one term in the expression with no `ga` and no `gascon` in it. Gravity
cancels everywhere except there, since `dqvi` already carries 1/g, so it should
scale with the atmospheric scale height, which here is 0.766 of Earth's. At a
representative LWP of 100 g/m2 that takes `ztau` from 30.8 to 22.6: **a 27% cut in
cloud optical depth, in the direction of a darker planet.** `0.00021`, a 0.21
g/m3 reference in-cloud liquid density, has the same standing and no dependence
on anything.

`inherited-earth-constants.md` finding 2 (PHYS-11) is the right neighbour and does
not cover this. It names the namelist multipliers `tswr1/2/3`, `acllwr`, `rcl1`,
`rcl2`, `acl2` and `clgray`, and it is about how a cloud reflects. How much water
the cloud is given has never been looked at, and the Stephens (1978) fit that
consumes it is equally uncited.

Beside it, the stratiform trigger `rcrit(:)=MAX(0.85,MAX(sigma(:),1.-sigma(:)))`
(`rainmod.f90:92`) is a bare floor and a bare profile function. At the eight
interior levels `1/(1-rcrit)^2` = 44, so a 0.01 move in `rcrit` moves cloud
fraction by 0.01 to 0.03 absolute wherever the modelled relative humidity sits
near threshold. ExoPlaSim ships `rcritmod` and `rcritslope` for exactly this and
neither appears anywhere outside `vendor/`.

## 9. The boundary-layer scheme runs entirely at compiled defaults

Every run's `fluxmod_namelist`, `seamod_namelist` and `surfmod_namelist` are
empty -- `&fluxmod_nl /END` -- so `fluxmod.f90:25-30` is live and unexamined:

    zumin      = 1.    ! minimum wind speed for PBL exhcange (m/s)
    vdiff_lamm = 160.  ! const. used in vdiff (see parameterization)
    vdiff_b    = 5.    !        "
    vdiff_c    = 5.    !        "
    vdiff_d    = 5.    !        "

`vdiff_lamm` is the asymptotic mixing length in metres, the single number setting
free-tropospheric vertical diffusivity, and `zlamh = vdiff_lamm*SQRT(3*vdiff_d*
0.5)` at `:783` derives the heat version from it. The three fives are the Louis
stability-function coefficients, cited at `:243` to "ECHAM REPORT 218" -- a
citation to a report, not a derivation, and none of the three is separately
justified. `zumin` is documented as m/s and applied at `:204` as a floor on
`u^2+v^2` and again at `:816` as a floor on a wind difference: one constant, two
quantities, one wrong unit in the comment. `parameter(ztscal=250.)` at `:762`, the
reference temperature converting sigma to metres for the mixing-length profile,
has no comment at all.

`fluxmod.f90:254-255` carries a hidden `g^(1/3)`:

    zrifh(jhor)=(1.+(0.0016*zdth**(1./3.)/SQRT(zabsu2(jhor))/zkblnz2)**1.25)**0.8

The Miller et al. (1992) free-convection enhancement, on the `zri <= 0 .and. dls
< 1` branch, which is unstable ocean and therefore most of the ocean most of the
time. In the strongly convective limit the bracket collapses to
`0.0016*dth^(1/3)/(|U|*C_N)`, which is the convective velocity scale with `g` and
the inversion height folded into the coefficient, and there is no `ga` in the
expression. At 12.81 against 9.81 the coefficient is low by `1.306^(1/3)` = 1.093,
so **free-convection latent and sensible exchange over calm unstable ocean is
about 9% weak**, over the warm ocean where this world's evaporation is.
`missed-couplings.md` reaches `dtransh` and stops before this branch.

## 10. Earth's radius in `oceanmod`, under a scheduled measurement

`oceanmod.f90:24`, used at `:1344` as `zfac=hdiffk(jlev)/plarad/plarad`:

    parameter(PLARAD=6.371E6)         ! Earth radius (m)

`oceanmod` uses `resmod` rather than `pumamod`, so this is its own parameter and
shadows nothing. It is dormant: `NHDIFF=0` in every run namelist. It is recorded
anyway because of what `config/planet.yaml` says about that key -- the horizontal
diffusion keys "exist to be BRACKETED in an A/B, not tuned into the baseline",
which is `clim-65`. When that A/B runs, `hdiffo` divides the requested `hdiffk` by
Earth's radius squared on a planet of 1.20 radii, making the effective diffusivity
1.44 times whatever the arm declares, with nothing in the output to say so. The
bracket would come back wrong. This is a trap under a scheduled measurement rather
than a dormant defect.

`ocean-and-marine-biosphere.md:95` catches the identical `parameter(radea=
6.371E6)` in `cpl.f90`, which is genuinely dead; the live twin two files away was
missed.

## 11. A fourth copy of Earth's lapse rate, in the postprocessor

`pyburn.py:217`:

    RLAPSE      = 0.0065    #International Standard Atmosphere temperature lapse rate in K/m

Used at `:2001, :2004, :2841, :2844` to build sea-level pressure, code 151 /
`psl`, which is in both `REGULAR_CODES` and `SNAPSHOT_CODES`
(`run_exoplasim.py:31,41`).

`inherited-earth-constants.md` finding 3 (PHYS-12) swept three copies of 6.5 K/km
and built `lib/lapse.py` to measure the rate from the model's own profile, 7.8
K/km per `clim-34`. **This is a fourth copy the sweep did not reach, and it is the
only one that writes a field into every climatology.** 6.5 against 7.8 is 20% low
and the reduction error grows with terrain height. Nothing downstream currently
reads `psl`, which is the only reason this is not already wrong in a quoted
number. `alpha = gascon*RLAPSE/gravity` does use the correct gravity, so the lapse
rate is the entire defect.

Beside it, `pyburn.py:2003` and `:2843`:

    tstar[tstar<255.0] = 0.5*(255+tstar[tstar<255.0])

ECMWF's cold-surface guard for the reduction, calibrated on Earth's surface
temperature distribution. It fires wherever the modelled bottom-level temperature
is below 255 K, which on a 32-degree-obliquity world is a large fraction of the
land area for much of the orbit, and it warms `tstar` by up to `(255-T)/2`, which
is 10 K at 235 K. Neither the threshold nor the halving is derived and there is no
config key.

## 12. The vertical grid comes from a library subclass default

`run_exoplasim.py` instantiates `exo.Earthlike` and passes neither `vtype` nor
`modeltop`. `Earthlike.configure` (`__init__.py:4409`) supplies `vtype=4` and
`modeltop=50.0`, so `NEQSIG=4` and `PTOP=5000.0` Pa reach the model from a library
subclass default. Neither has a key in `config/planet.yaml`.

`parameter-decisions.md:558` records "`modeltop` None, giving `PTOP` 5000 Pa,
appropriate for a 1 bar atmosphere". **The premise is false.** If `modeltop` were
None then `if modeltop:` is false, `PTOP` is never written, and the model uses
`plasimmod.f90:315`'s compiled `ptop = 7500.0` Pa. The number is right only
because `Earthlike` supplies 50 hPa. The same document opens "configure() takes
102 parameters. We set 26; the rest sit at defaults" -- two of those defaults are
the subclass's, and `vtype` is not mentioned in that audit at all.

`vtype` is not cosmetic. `plasim.f90:1716-1725` (`neqsig==4`) rescales the sigma
polynomial so the top half-level sits at `ptop/psurf`; the `neqsig==0` fallback at
`:1752` runs the same polynomial unrescaled, giving sigma 0.0766. The two put the
model top at 5000 Pa against 7660 Pa.

The polynomial itself, `sigmah = 0.75*zsk + 1.75*zsk**3 - 1.5*zsk**4`
(`plasim.f90:1720`), places all ten levels and has three coefficients with no
comment, no citation and no derivation anywhere in the tree. The only property
recoverable from the code is that its derivative vanishes at the surface.

## 13. `gamma` is 0.01 at T21, and the note records the T42 branch as live

`rainmod.f90:31` declares `gamma = 0.01`; `:88-90` sets `gamma = 0.007` only when
`NTRU==42 .and. NLEV==10`. `config/planet.yaml` reads `resolution: T21`, `layers:
10`, so the branch does not fire. `resolution-tuned-parameters.md:88` states the
T42 branch "which this project's configuration does hit". The documentation no
longer matches the code.

`gamma` is the fraction of the sub-saturation deficit that falling precipitation
evaporates per timestep (`rainmod.f90:2244-2246`; it is divided by `deltsec2`, so
it is not a rate and carries no time unit). It controls how much precipitation
re-evaporates before reaching the ground, which is P-E over land directly. The
two branches differ by 43%.

This is not a one-off. The resolution has flipped between T21 and T42 five times
and `resolution_timestep_minutes` gives both rungs the same 45-minute step, so
nothing compensated. A precipitation-physics constant has been oscillating under a
key that nothing connects to precipitation. Neither branch has a derivation
upstream either.

## 14. `OMP_STACKSIZE` is in every instrument and in none of the model

Nine gate and bench scripts set `OMP_STACKSIZE=512M` with `ulimit -s unlimited`:
`verify_shtns_model.sh:140`, `verify_threaded_numerics.sh:137`,
`verify_shared_determinism.sh:106`, `verify_weight_factorisation_model.sh:103`,
`verify_omp_collectives.sh:63`, `profile_transforms.sh:48`, `profile_memory.sh:64`
and `:70`, `attribute_barrier_wait.sh:41`, `bench_ab.py:62`. `bench_ab.py:74-79`
records that without the `ulimit`, at T127 the master overruns a 16 MB limit and
segfaults.

The production launcher sets neither. `__init__.py:522-523` builds `_exec` as
`OMP_NUM_THREADS=%d OMP_PLACES=cores OMP_PROC_BIND=close ./`, and a grep for
`OMP_STACKSIZE` or `ulimit` across `__init__.py`, `run_exoplasim.py` and
`continue_exoplasim.py` returns nothing. **The one setting that decides whether a
threaded T127 or T170 run survives at all exists only in the thing that measures
the model and not in the thing that runs it**, and the threading work is headed
at those rungs. 512 is also a round number: `threads-instead-of-ranks.md:61` says
only that the largest spilled local is 14 MB, so it must be "well above the 8 MB
default".

## 15. T31 is unbuildable, and the rung table is spelled four more times

`CMakeLists.txt:57` hardcodes `require_one_of(PLASIM_NLAT ... 32 64 96 128 160 192
256)`, and the help string at `:31` repeats it. `lib/rungs.py` has `T31: 48`.
`build_model.py --res T31` resolves cleanly through `rungs.RUNGS` and then dies at
`cmake` with `FATAL_ERROR`. It fails loudly rather than silently, but T31 is
unbuildable and nothing says so.

This is the defect `lib/rungs.py`'s own docstring says it was created to end
(SPAT-2: "four scripts each carried their own copy of the mapping, and two of the
four were missing rungs the others had"). Four more private copies survive:

| file | what it carries |
| --- | --- |
| `cache_budget.py:32-40` | omits T63 and T106, invents **T213** (213, 320, 640), which is on no ladder anywhere; its comment claims "as `build_model.py` resolves them" where `build_model.py:58` is `RESOLUTIONS = dict(rungs.RUNGS)` |
| `verify_transform_roundtrip.sh:30-31` | missing T63 and T106; its two sibling scripts already shell out to `rungs.geometry()` |
| `stability_probe.py:63` | five of eight rungs |
| `filter_timestep_matrix.py:103-104` | cites SPAT-2 in the comment directly above its own local copy |

`cache_budget.py` uses `choices=sorted(RESOLUTIONS)`, so it refuses T63 and prints
a per-die budget for a rung with no grid and no binary. It is the instrument for
the 32 MB working-set rule.

## 16. The albedo bracket can no longer bracket

`compare_albedo_bracket.py:119`:

    name = "lithology" if mean > 0.256 else "vegetated"

The docstring cites "bare rock at 0.315 land-mean albedo and vegetated at 0.223".
Measured 2026-08-24: `inputs/t42/albedo_report.json` gives
`land_mean_bare_rock = 0.2483`, below the classifier threshold. Every `--mode
lithology` run on the current build therefore classifies as `"vegetated"`, and
since `temps` is a dict keyed on that name (`:163`) the two arms of a bracket
collapse to one entry and the script reports "cannot bracket". This decides the
verdict on the experiment that establishes whether the vegetation feedback is
first-order.

## 17. The energy fixer's controller state is not in the restart

`denergyfix`, `denergyacc`, `nenergyacc` and `nenergywin` appear in no
`put_restart_*` call (`plasim.f90:902-1012`) and no `get_restart_*` call
(`:1163-1260`). Every segment therefore begins with `denergyfix = 0` and also
discards its first window, so the first `2*ntspd` = 80 steps of every segment run
with no correction and then step to the full value.

At 5850 steps per orbit that is 1.4% of a segment, worth roughly 0.011 W/m2 on the
segment mean against a 0.05 W/m2 convergence criterion: small, one-signed, and it
makes each segment's first two model days non-comparable to the rest. The sharper
consequence is that **any run shorter than 80 steps never applies a correction at
all**, which includes the 20-step SHTns verification gate. Nothing in that gate
exercises the fixer.

Adjacent to `world-5qy`, which is about `naccuout` and the output boundary, but a
different object.

## 18. The diffuse aerosol two-stream is the direct one

`radmod.f90:2661-2664` for band 1, `:2693-2696` for band 2:

    zaertf1s(:,jlev) = MIN(25.,(ztemp1(:)*aod1(:,jlev))/zmu00)  ! ... using zmu00 not zmu0!
    zaerd1s(:,jlev) = (((zaeru1(:)+1.0)**2.0)*EXP(zaertf1(:,jlev)) - ...
    zaert1s(:,jlev) = (4.0*zaeru1(:))/zaerd1(:,jlev)

`zaertf1s` and `zaerd1s` are assigned and never read: their only occurrences in the
file are the declarations at `:2354` and `:2356` and the assignments at
`:2661-2662` and `:2679-2680`. `zaerd1s` is built from `zaertf1`, the direct-beam
path, and `zaert1s` and `zaerr1s` then divide by `zaerd1`, the direct-beam
denominator. So the diffuse-beam aerosol transmission and reflection are
bit-identical to the direct-beam ones, the comment describes something that does
not happen, and the factor-of-two diffusivity the scattered stream should carry is
absent.

Inherited from upstream, but the fork rewrote these exact lines for the per-cell
arrays and preserved it, then wrote the new conservative-scattering branch to set
`zaert1s = zaert1` and `zaerr1s = zaerr1` explicitly at `:2681-2682`, which makes
the defect now read as deliberate. Latent: `NDUSTRAD=0` and no run has ever had
aerosol radiation on. It is the path `clim-40` and `dust-13` switch on.

## 19. A second radiative feedback, 11% off the canonical slope

`assess_convergence.py:69`:

    FEEDBACK_W_M2_K = 1.31                            # measured, not assumed

with no pointer to where. Measured 2026-08-24: `lib/sensitivity.py`'s
`SLOPE_K_PER_FLUX_RATIO = 202.0` gives `kelvin_per_w_m2(alpha=0.3)` = 0.84812, an
implied radiative damping of **1.1791 W/m2/K** -- the 1.18 that
`forcing-bundle-predictions.md:370-372` uses and attributes to that module.

CLAUDE.md names `lib/sensitivity.py` as the one flux-to-kelvin conversion.
`assess_convergence.py` does not import it. `FEEDBACK_W_M2_K` sets `tau_expected`,
which sets the fallback remaining-offset judged against `OFFSET_TOLERANCE_K =
0.15`, so it is an input to a pass/fail verdict. Beside it at `:68`,
`SLAB_HEAT_CAPACITY = _MLD * 1025.0 * 3990.0` is uncited where the model's own
values, cited to source in `close_state_energy.py:67-68` and
`close_ocean_energy.py:121-122`, are `CRHOS = 1030.0` and `CPS = 4180.0`, 5.3%
higher. The slab is modelling the model's mixed layer, so the model's constants
are the right ones.

## 20. The transcription route, and three gates that cannot fail

Two structural patterns, each with several instances.

**Config values are hand-transcribed from artifacts and nothing compares them.**
The generators print a value and a human retypes it into `config/planet.yaml`.
`check_consistency.py` never opens `shortwave_band_weights.json` or
`cloud_band_weight.json`, and `lib/provenance.py` only lists the keys as inert or
active. The values agree today -- `h2o_sw_weight` 1.346, `co2_sw_weight` 1.510,
`cloud_absorption_scale` 1.192, each verified against its JSON on 2026-08-24 --
but the route has already drifted twice: finding 5 in the model source, and the
dust tables, where `dust_forcing.py:139-141` and `dust_optics.py:88-91` now read
`vegetation_albedo` 0.165 and `playa_clastic` 0.23 from config while
`notes/dust.md:224-227, :325-328, :396-400` and `analysis/dust_optics.json` still
say 0.18 and 0.40. Interpolating that note's own near-linear table, playa TOA net
falls from +4.3 to roughly +0.8 W/m2 and the land half of the global mean from
about +2.05 to +0.15, roughly 13 times, against the 1.5 W/m2 reopening threshold
at `dust_forcing.py:472`. `dust_optics.json`'s `"warms_over": ["playa fill", ...]`
becomes false at 0.23, its `critical_surface_albedo` being 0.3001.

**Three gates have bounds several decades from anything ever observed.**
`verify_threaded_numerics.sh:70-71`, copied verbatim to
`verify_shtns_model.sh:59-60`, sets `BIRTH=1e-11` and `JUMP=1e4` as hard `rc=1`
gates; the 50-line docstring above them derives the run lengths carefully and
neither bound. The observed jump in the only run on record is 5.4x, four decades
below the gate, which matters because `shtns-viability.md:907` leans on "the jump
bound already tests the growth's shape". `transform_exactness.py:57` sets
`ROUNDOFF = 1e-10` where double-precision round-off on a global mean over 64x128
is of order 1e-14, and it decides whether the `world-bxr` energy-sink measurement
is valid or void; its neighbour `STORAGE = 1e-6` is properly derived from float32
output. `corrk_cross_check.py:206-234`'s `run_checks` contains no threshold, no
comparison and no raise, while its docstring claims it proves the division works
to 0.2% and `--checks` returns at `:309` before any number is computed.

---

## Smaller findings, recorded so they are not re-derived

| where | constant | what is missing |
| --- | --- | --- |
| `radmod.f90:2504` | `35./SQRT(1.+1224.*mu^2)` | Kasten airmass; `1224 = 35^2-1` and 35 is Earth's horizon airmass. Recomputed for this world: 43.1 and 1854. Effect is small (1.014x at mu = 0.1) and the arithmetic is here so it is not redone |
| `radmod.f90:2266-2268` | `aa=0.2542857142857143`, `bb=0.8229693877551021`, `c0=0.14997959183673468` | three sixteen-digit literals, no comment, no source; dead at `newrsc=0` and live on one namelist key |
| `radmod.f90:2553, 2554, 2761, 3201` | `101100.0` | `psurf` retyped as a literal in four places; this project sets `PSURF = 99999.99999999999`. `:2761` is live and understates surface Rayleigh reflectance by 1.1% relative, about 0.11 W/m2 |
| `radmod.f90:81-82, 1987-1993` | `bo3=20000.`, `co3=5000.` | ozone layer placed at a fixed geometric altitude; `ozone.md:15` names both and states neither their unit nor their planet dependence. 20 km is p/p0 = 0.093 on Earth and 0.045 here |
| `radmod.f90:199, 449-450, 576` | `minwavel = 316.036116751` | a twelve-digit value whose only stated justification is that it reproduces the hardcoded `zsolar1 = 0.517`; `:576` then uses `zsolar1` as the Sun's reference 67 lines before the same routine overwrites it with this star's |
| `radmod.f90:264` | `dusthsc(NAERSP) = 3000.0` | fork-introduced as an array; every species 2 to 4 silently inherits species 1's scale height, with no guard, where `dustqlw` aborts at zero. Duplicates `aeolian/config/dust.yaml`'s declared 3000.0 |
| `seamod.f90:154-158, 269-272` | `0.025` per K, anchored at `273.` | the sea-ice albedo ramp. Its endpoints ARE star-corrected (`radmod.f90:725-753`); the slope between them was not, and the anchor is neither `TMELT` (273.16) nor `TFREEZE` (271.25) |
| `landmod.f90:406-410, 557-561` | `0.01` | snow-cover masking depth, 1 cm water equivalent, with no dependence on the surface being covered, on a world carrying a lithology roughness field spanning 0.025 to 11.2 m. Not covered by `grav-8`, which is about density |
| `landmod.f90:54, 98` | `dztop = 0.20`, `dsoilz = 0.4,0.8,1.6,3.2,6.4` | no derivation. The diurnal damping depth here is 0.161 m against Earth's 0.144, so `dztop` is 1.24 damping depths against Earth's 1.39 and the modelled land diurnal amplitude is systematically larger |
| `rainmod.f90:29` | `pdeepth = 70000.` | Earth's 700 hPa. The sigma is right by accident on a 1 bar surface; the height is 2.7 km here against the 3.5 km the threshold was chosen for |
| `rainmod.f90:1905-1906, 1926-1928` | `zcca=0.245`, `zccb=0.125`, `zccmax=0.8`, `zccmin=0.05` | convective cloud fraction from a log of a dimensional rain rate. `solar_day` IS correctly planet-aware, so the conversion is right; what is missing is any record of which day convention the Earth fit was made under, which leaves the sign of the correction undecided |
| `landmod.f90:52` | `drhsfull = 0.4` | `missed-couplings.md:194` identifies it as the entire land-evaporation channel at `NVEG=0` and never asks where 0.4 came from |
| `glaciermod.f90:47` | `rhoglac = 850.` | firn rather than glacial ice (917), and nothing says which was intended. A third density in the chain `grav-8` lists two of |
| `landmod.f90:1185`, `glaciermod.f90:303` | `zcvel=4.2`, `zcexp=0.18` | river routing, duplicated in two files that both run; the slope is taken in geopotential, so the hidden `g` dependence is `1.306^0.18` = 1.049 |
| `gaussmod.f90:46, 48` | `NITER=100`, `ZEPS=1.0e-15` | no post-loop convergence test and no diagnostic, so a node that fails to converge produces a plausible weight silently. `1.0e-15` is also a default-real literal in a `real(kind=8), parameter` |
| `plasim.f90:1437-1439, 1460-1463` | `tfrc` = 20, 30, 100 days | Rayleigh drag on the top layers, live at `NLEV==10`, never in a namelist, never derived, and by finding 2 applied as 25 and 125 planetary days |
| `plasim.f90:3689` | `abs(zfixc) > 2.0` | the energy fixer's per-window clamp; fork-introduced, states its role and no criterion, where the divergence bound eleven lines later does state one |
| `plasim.f90:2219-2226`, `plasimmod.f90:311-313` | `alr=0.0065`, `tgr=288.0`, `dtrop=12000.0` | the cold-start initial atmosphere, all Earth's, all echoed live. Bounded because a bootstrap equilibrates away from its initial state, so this costs spin-up length and blow-up risk rather than equilibrium. `clim-71` proposes to cold-start through it |
| `calmod.f90:21-22` | `m_days_per_year=360` | behind a comment saying the values are copied from pumamod in `calini`, which does not copy them. Gives a calendar year of 5760 steps against an orbit of 5850. Reaches output metadata only: the orbital forcing takes its phase from `mod(nstep,n_steps_per_year)` at `radmod.f90:1848` |
| `plasimmod.f90:452` | `landhoskn0 = 15.0` | carries its own admission that it is tuned to T21, with no branch for any other rung. Unreachable at `physics_filter: "gp|exp|sp"` |
| `__init__.py:2063` | `oceanzenith="ECHAM-3"` | sets `NECHAM=1`, the switch `ocn-7` and `ocn-22` describe from the Fortran side. There is no `NECHAM` string anywhere under `exoplasim/`: it reaches the model from a Python keyword default |
| `__init__.py:2653` | `NSTPW` = 7200 min | five Earth days; 5850 % 160 = 90, so 90 timesteps accumulate unwritten under `--low-io`. `run_stellar_cycle.py:212-217` works around it and `run_exoplasim.py` does not |
| `__init__.py:2711` | `radius*6371220.0` | ExoPlaSim's Earth radius, where Orogen declares `radiusKm: 7645.2` (1.20 x 6371.0) and this project's own scripts use 6371 km. 264 m of radius, 7e-5 in area: negligible physically, a rule-2 shape structurally |
| `pyfft.f90:149`, `pyfft991.f90:148` | `exp(-8*(n/NTRU)**8)` | the postprocessor's own physics filter, carrying the superseded gamma of finding 4. Inert at `physfilter=False`, live the moment anyone enables it |
| `__init__.py:2754` | `outgassing = 50.0` | written as `VOLCANCO2` into every run's `carbonmod_namelist` regardless of `NCARBON`, so every run directory records an outgassing rate nobody chose. `volc-8` is open on this world's outgassing requirement |
| `dust_forcing.py:283, 616` | `lambda x: np.ones_like(x)` | a flat spectral weight where `main()` passes the star. Band-2 mass extinction efficiency is 0.5759 flat against 0.7646 weighted, -24.7%, and band 2 carries 61.6% of this star's flux |
| `dust_forcing.py:100, 107-108` | `TRANSMISSION=0.79` (applied squared), `WINDOW_TRANSMITTANCE=0.7` | a 0.624 multiplier on every shortwave forcing number in the component, and a declared-but-never-computed bracket on the longwave half |
| `shortwave_band_weights.py:204` | `EARTH_CO2_SHORTWAVE_W_M2 = (1.5, 2.5)` | called "THE CHECK THAT CAN FAIL" in the docstring, repeated in three notes as what the literature measures, with no source and absent from `references/INDEX.md`. The range is 67% wide |
| `shortwave_band_weights.py:695` | `0.65/0.35` above/below cloud | the bracket ends are argued and the collapse is not; `shortwave-water-vapour.md:255` admits it and the code does not. The same 0.65 is applied to CO2, where `shortwave-co2.md:190-196` says outright it is wrong |
| `cloud_band_weight.py:81` | `SIGMA_G = 1.35` | bare, in the line below `R_EFF_UM` declared as "Bracketed, not chosen". It fixes the modal radius the Mie integral runs at, and the reported bracket spans only `r_eff` |
| `run_exoplasim.py:27`, `analyze_climatology.py:30`, `run_stellar_cycle.py:54`, `shortwave_band_weights.py:715` | 365.2568983 / 365.2425 / 365.25 / 86400*365.25 | four spellings of the year outside `lib/orbit.py`, whose docstring says duplicated year literals "drift the moment the flux changes, which is exactly what happened". `shortwave-water-vapour.md:230`'s 974 mm/yr is therefore per Earth year on a 182.8-day world: `CLIM-25` in a different file |
| `bench/bed_check/bed_manifest.json` | `"steps": 60` | against `make_profile_bed.py:152`'s cited default of 600. At T42's 18 ms/step that is 1.1 s of integration against a 1.6 s fixed startup, so every bucket share this bed reports is startup-weighted. `failure-modes.md` class 34 |
| `bench_ab.py:145`, `bench_rank_layout.py:296`, `sweep_compiler_flags.py:129`, `shtns_variant_sweep.py:177` | 6 / 6 / 8 / 4 rounds | every one of these files documents interleaving and rotation rigorously and none says how many rounds are enough. `symmetric-transforms.md:152` records the count being changed after seeing a result |
| `compare_restarts.py:66` | `--abs-floor 1e-12` | can set `r = 0.0`, turning a DIFFERENT into an identical. The nine-line help argues the purpose and never the magnitude, and its claim that model quantities are O(1) to O(1e5) is contradicted at `:122` in the same file, which records `dql` at 1e-9 |
| `verify_shtns_equivalence.sh:66` | `-O2 ... -finit-real=zero` | a hand-written flag line for a probe compiling the same `legmod`/`shtnsmod`/`gaussmod`/`fftmod` the model does, without `-march=znver4` (so no FMA contraction, i.e. different rounding in the reduction loops under test) and with the very flag `the-zeroing-is-an-init-flag.md` deliberately removed for a 26% gain. `build_model.py:160-165` states that `--extra-flag` is the only undeclared route to the compiler |
| `CMakeLists.txt:101` | `set(CMAKE_C_FLAGS "-O3")` | the only flag CMake supplies itself, against the file's own header at `:21-22` saying the flag line arrives whole from config and CMake must not add to it. It also skips the `-ffile-prefix-map` pair, so C objects sit outside the byte-comparability argument |
| `sweep_compiler_flags.py:57-71` | `"drop": ["-fcheck=all", ...]` | six of seven arms drop a flag the `production` profile does not carry, so `build_model.py:151-155` raises and the sweep catches it: the flag sweep currently reduces to a single `stock` arm |
| `build_model.py:58-59` | `NEEDS_FFT991 = {"T63","T106"}` | the set is right and the stated reason is not: T127's NLON is 384 and T31's is 96, neither a power of two, and both go through `fftmod`. The real criterion is membership of the two `nallowed` tables |
| `build_surface_soil_water.py:188` | `model.get("lake_dwmax_m", 0.2)` | `lake_dwmax_m` is not in `config/planet.yaml` at all, so the fallback is what runs, and the `.get()` reads as though config supplies it |
| `build_surface_albedo.py:163-167` | `--tree-albedo 0.13`, `--grass-albedo 0.19` | two of the three constants `inherited-earth-constants.md` finding 1 named in one sentence; the third was re-weighted to this star and these were not |
| `run_exoplasim.py:1885-1887`, `run_stellar_cycle.py:189-191` | `fixedorbit=True, keplerian=True, meananomaly0=0.0` | hardcoded in two drivers, in no config block, with zero hits in `docs/`, `notes/` or `bd`. With `longitude_vernal_equinox_degrees` they set the seasonal phase that `CLIM-23` covers only half of |
| `p_earth.f90:41` | `akap = 0.286` | live, feeding `acpd = gascon/akap`, where `GASCON` was derived from the declared composition and this was not. Correct value 0.28562, implying cp 1004.90 against 1003.56: 0.13%, worth 0.017 K/km on the dry adiabat |

## Checked and clean

Recorded so it is not re-audited.

**No value from `config/planet.yaml` leaks into the Fortran.** Gravity, radius,
rotation, flux, effective temperature and CO2 were each searched for across
`vendor/exoplasim`, `exoplasim/`, `config/` and `lib/`; the only hits are the
config file itself and the registry's prose. The full config-to-namelist chain
traces end to end for `gravity_m_s2`, `radius_earth`, `rotation_hours`,
`baseline_flux_earth`, the four partial pressures, `eccentricity`,
`obliquity_degrees`, `longitude_vernal_equinox_degrees`, the derived orbital
period, `mixed_layer_depth_m`, `filter_kappa`, `filter_power`, `salinity_psu`,
`ozone_scale`, the five shortwave weights, the two trace-gas abundances and
`cold_start_seed`. `N_DAYS_PER_YEAR` and `N_RUN_STEPS` are both correctly
overridden past `configure()`'s Earth-360-day fallbacks.

**`p_earth.f90`'s remaining constants are either overridden or inert.** `alr =
0.0065`, `tropical_year = 31556956.0` and `nplanet = 3` are set, MPI-broadcast and
never consumed anywhere in the live source. `solar_day` and `sidereal_day` are
recomputed from `rotspd` at `plasim.f90:1490-1503`, so Earth's 86400 does not
survive into a run and the convective cloud fit at `rainmod.f90:1926` gets the
right day length. `ra1`, `ra2`, `ra4` and `tmelt` are properties of water.

**The fork's spectral work cites everything it carries.** `SHTROOT`, `SHTRINV`,
`SHT_ORTHONORMAL`, `SHT_QUICK_INIT`, `zeps=0.0` and `shtns_use_threads(1)` are all
measured and cited in-file to `probe_shtns_*`, `shtns-viability.md` and
`shtns-algorithm-selection.md`. `mpimod_omp.f90`'s only bound is derived.
`legmod.f90` after the two-matrix rewrite contains no magic numbers.

**The energy fixer's own constants carry their criteria in-code**: the `ntspd`
window, the discarded first window, the tendency-not-increment choice and the 100
W/m2 divergence bound are all argued at `plasim.f90:3661-3703`. Only the 2.0 W/m2
clamp is not.

**Gravity is handled correctly in three places it could have been missed**:
`radmod.f90:1928`'s explicit `(9.80665/ga)` in the Rayleigh column mass,
`pyburn.py`'s `alpha = gascon*RLAPSE/gravity`, and `seamod.f90:262-263`'s Charnock
term, which computes `charnock*|tau|*gascon*T/(ga*dp)` with the model's own `ga`.
`seamod.f90:19`'s comment on `charnock` is copy-pasted from `albsea` four lines up
and describes the wrong quantity, but the code is right.

**One diagnostic that is wrong and reads as authoritative.** `print_planet` is
called from `plasim.f90:196`, before `readnl` recomputes `sidereal_day` from
`rotspd` at `:1490`. Every run log's planet table therefore reports Earth's mass,
volume, equatorial and polar radii, ellipticity, density, Bond albedo, black-body
temperature, perihelion and aphelion under the header "Simulating: Earth", with a
correct mean radius one row below Earth's equatorial one, and two derived rows --
a 23.9345 h rotation and a 183.3019-day orbit -- computed before the rotation is
applied. The model then integrates a 30 h rotation and a 146-day orbit. Nothing
in the physics reads any of it.

## What this audit did not cover

Stated so the coverage is not overread. Dormant modules were checked for liveness
and not mined: `simba.f90`, `carbonmod.f90`, `tracermod.f90`, `trc_routines.f90`,
`newsnow.f90`, `buildice.f90`, `hurricanemod.f90` and the LSG ocean, all of which
`dormant-exoplasim-modules.md` already classifies. The deselected convection
schemes `rainmod_bm.f90`, `rainmod_mca.f90` and `rainmod_kuo_old.f90` carry
identical copies of findings 8 and the `rcrit` note and were not separately
listed. I did not audit the Lacis-Hansen and Sasamori coefficient tables beyond
confirming they are covered by `ozone.md`, `shortwave-water-vapour.md` and
`corrk-cross-check.md`; the LPJ-GUESS, hydrography, pedology or minerals
components; or Orogen's generation constants.

## Tasks

Tracked in the `bd` issue tracker under the `opaque-const` label, not restated
here.

| finding | id |
| --- | --- |
| 1. hyperdiffusion reaches one model level of ten | `world-720` |
| 2. every damping timescale is applied 1/rotspd too long | `world-rt1` |
| 3. `t0 = 250.0` and the energy sink that scales with it | `world-bmf` |
| 4. the transform gate certifies a superseded filter | `world-tez` |
| 5. the CO2 shortwave fit diverges from its patch and artifact | `world-vej` |
| 6. surface longwave emissivity | `world-qvu` |
| 7. the roughness anchor is computed at Earth's gravity | `world-44l` |
| 8. cloud liquid water and cloud optical depth | `world-ofn` |
| 9. the boundary layer runs at compiled defaults | `world-e2k` |
| 10. Earth's radius in `oceanmod`, under `clim-65` | `world-mll` |
| 11. a fourth lapse rate, and a 255 K clamp, in `pyburn` | `world-ld1` |
| 12. the vertical grid comes from a subclass default | `world-a05` |
| 13. `gamma` is 0.01 and the note says otherwise | `world-j6v` |
| 14. `OMP_STACKSIZE` is absent from the run path | `world-3nk` |
| 15. T31 is unbuildable and the ladder is spelled four more times | `world-b3s` |
| 16. the albedo bracket collapses to one arm | `world-ue5` |
| 17. the energy fixer's state is not in the restart | `world-fsr` |
| 18. the diffuse aerosol two-stream is the direct one | `world-sv7` |
| 19. a second radiative feedback, 11% off the canonical slope | `world-1n3` |
| 20a. the transcription route nothing compares | `world-2bf` |
| 20b. three gates that cannot fail | `world-134` |
| the solar/sidereal timestep, latent below dt 12.4 min | `world-r8o` |
| the run log's planet table | `world-1o4` |
| `clim-68`'s stale citations and inverted sub-claim | `world-ylw` |
| the compiled planet module, `akap`, and three inerts | `world-cwu` |
| `lake_dwmax_m` is not in config | `world-vus` |
| two vegetation endmembers left at Earth-Sun values | `world-9m5` |
