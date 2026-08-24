# Where the vendored model still describes Earth

Audited 2026-08-24, on `source_build: precarve-craton-10m` at T21 with the
staged namelists and run output of `run_2b20e3324bb0` and `run_2c42c68fe9ca`.
The subject is the SIMULATION of Vesper: a modelled atmosphere, a modelled
ocean, a modelled snowpack. Line citations are to
`vendor/exoplasim/exoplasim/plasim/src` unless another path is given.

`inherited-earth-constants.md` asked where an Earth number reaches this world
with nothing naming it, and hunted in the radiation's surface and cloud terms,
the orbit block and the flux derivation -- that is, in this project's own
configuration and analysis layer. `absent-and-inherited-physics.md` hunted the
ocean, the cryosphere and the soil parameters for the same thing. This one
takes the question into the FORK ITSELF: roughly 50,000 lines of Fortran and
15,000 of Python, swept by domain, with every load-bearing claim checked
against a staged artifact rather than against the source alone.

**This sweep opened twenty-nine findings: fourteen live and load-bearing,
fourteen live but bounded or diagnostic, and eleven dormant traps waiting on a
switch this project already intends to throw.** Each is kept below with its
evidence, and each now records what its fix left standing, because a finding is
the reason a change was made and outlives the change.

Two are open in full. Finding 14, one Earth-average soil thermal pair applied to
every lithology, whose six constants became namelist scalars but which still has
no per-cell field behind it (`world-sy9`). And the dormant `jtune` radiation
table, unreachable at every rung because `ndcycle` defaults to 1 (`world-ys9`).
Everything else below is either settled or reduced to a declared value at the
number the literal carried, which is stated finding by finding.

The individual constants matter less than the four mechanisms that keep
producing them, and those are stated first because they are what a reader needs
in order to predict where the next one will be.

---

## 0. A converged climatology was integrated on the wrong diffusion operator

Not an Earth-centrism finding. Recorded here because this sweep is what found
it and because it invalidates an artifact the tree currently trusts.

`continue_exoplasim.py` imported and called five of the six `declare_*` helpers
after `configure()` and omitted `declare_hyperdiffusion`. `configure()` rewrites
`plasim_namelist` wholesale on every continuation, so each continued segment
dropped back to PlaSim's compiled Earth tuning. The same two drivers also
omitted `filterpower`, so `NFILTEREXP` reverted from the derived 16 to the
library default 8.

**Closed by world-8bs.** `continue_exoplasim.py` now passes `filterpower` into
`configure()` and calls `declare_hyperdiffusion` beside the other five, so a
continuation restates the filter and writes the hyperdiffusion keys. The two
classes are handled differently on purpose and the comment at
`run_exoplasim.py` says why: `configure()` writes `FILTERKAPPA` and
`NFILTEREXP` unconditionally from its own defaults, so a dropped key there is a
WRONG value, while it does not write `NDEL`, `NHDIFF` or `TDISS*` at all, so a
dropped key there is an ABSENT one.

Checked across all nine run directories on 2026-08-24. One is missing every
key:

    run_2b20e3324bb0   NHDIFF=-    NDEL=-   TDISSZ=-        NFILTEREXP=8
    run_2c42c68fe9ca   NHDIFF=16   NDEL=4   TDISSZ=1.1143   NFILTEREXP=16
    run_52349d4bfb75   NHDIFF=8    NDEL=4   TDISSZ=2.2285   NFILTEREXP=16
    ... six more, all complete ...

`exoplasim/runs/INDEX.json` records that run as `climatology_complete` and
`converged: true`. At T21 the compiled fallback is `nhdiff=15`, `ndel=2`,
`tdissd=0.20`, `tdissz=1.10`, `tdisst=5.60`, `tdissq=0.1`, because the
`NTRU==42` branch at `plasim.f90:1544-1551` does not fire at T21. So that run
integrated with a fourth-order operator instead of an eighth-order one, a
damping cutoff at 0.71 of the truncation instead of 0.38, and humidity damping
7.4 times too strong.

`expected_namelist_keys()` in `run_exoplasim.py` exists for exactly this class
and covered the shortwave gas weights, `O3SCALE`, `TFREEZE`, `NENERGY`,
`NENER3D`, `NCONVTIME`, `NDEALIAS`, `NENERGYFIX` and `PNU`. It did not cover
`NHDIFF`, `NDEL`, `TDISS*`, `NFILTEREXP` or `FILTERKAPPA`, so
`verify_staged_namelists` ran on every continuation and printed "namelists
verified" over both reversions. It covers all six now, along with the keys the
later findings below declared.

## The four mechanisms

**I. The build compiles Earth's planet module.** It is the ONLY planet module:
`p_earth.f90`. `PLASIM_PLANET` and the unbuildable `p_exo.f90` and `p_mars.f90`
beside it were removed under world-58v, and the three keys `p_exo.f90` alone
exposed -- `alv`, `als` and `tmelt`, the latent heats of vaporisation and
sublimation and the melting point -- are now in `p_earth.f90`'s `planet_nl` at
the values `plasimmod.f90` already declared, broadcast to every thread. What
remains Earth-shaped here is the module's identity, `yplanet="Earth"` and
`nplanet=3` with the Earth Jan-1 `meananom0`, which the run's own diag prints
as `YPLANET="Earth"`.

**II. `day_24hr` carried two roles at once, and one of them was the model's
unit of time.** `day_24hr = 86400.0` (`plasimmod.f90:958`) is hardcoded, appears
in no namelist group anywhere in `plasim/src`, and nothing in the source ever
assigns it. It is right for the role of converting a timescale ENTERED in
24-hour days into seconds, which is the unit `config/planet.yaml` derives its
`timescales_days` in and the unit `dayseccheck` decides against. It was also
being used as the model's unit of time, which is `1/ww = sidereal_day/TWOPI`,
and the two coincide only on Earth, where `rotspd` is 1. **world-rt1 separated
them by name.** `restim`, `tdissd`, `tdissz`, `tdisst`, `tdissq` and `dampsp`
are nondimensionalised with `sidereal_day` at `plasim.f90:1969-2002` and
`:1751-1753`; `deltsec`, `mtspd` and `dayseccheck` keep `day_24hr`, which is
the unit their values are entered in. The comment at `plasim.f90:1620-1633`
carries the distinction at the point where both are in scope.

**III. Deleting Earth's input files does not yield a Vesper default.**
`configure()` correctly removes the shipped Earth boundary fields, and the diag
confirms codes 169, 210, 211, 709 and 903 all print `Init ... internally`. The
model did not then fall back to something neutral: it CONSTRUCTED an Earth
Arctic ice field, or read an array nothing ever initialised. The absent-field
path now refuses or takes a declared value instead, finding by finding; the
mechanism is what to expect of the next field this project stops supplying.

**IV. Many of these constants had no namelist key at all.** A large share of
what follows sat in Fortran `parameter` statements or in module variables absent
from every namelist group. Those cannot be bracketed, swept, or declared in
`config/planet.yaml`, and moving one means a source edit and a rebuild of every
binary under rule 4. Most of the ones below are namelist keys now, at the value
the literal carried, which changes no result and makes each one an arm; where
the value itself is still Earth's the finding says so.

---

# Live and load-bearing

## 1. Every cold start begins from a constructed Earth Arctic ice field and an ocean at its freezing point

`iceini` left `xclsst` at its `-999` sentinel when code 169 was absent, and then
did `where (xclsst <= TFREEZE) xclicec = 1.0`, which puts full ice cover on
every ocean cell. `make_ice_thickness` (`icemod.f90:2079`) shapes it with
`hmaxn=3.0`, `hmaxs=0.50`, `cminn=0.1`, `cmins=0.25` and a fourteen-element
monthly modulation table credited in its own comment to "CCM3 report pp
127-129".

Separately, `oceanmod.f90` declared `real :: yclsst(NHOR,0:13)` with no
initialiser while its neighbours all had one, then did
`yclsst = MAX(yclsst,TFREEZE)` and copied it into the ocean temperature, so
whatever was in memory was clipped up to the freezing point.

**Closed by world-6fh, and the cold start is now declared.** `icemod_nl` carries
`tsst_eq`, `tsst_pol` and `hice_ini` (`icemod.f90:97` and neighbours), which
`config/planet.yaml`'s `ocean.cold_start` block sets: 300.0 K at the equator,
273.0 K at the pole, 0.0 m of ice. The default cold start is therefore ICE-FREE,
which is the side of the ice-albedo feedback a flux sweep could not otherwise
approach a branch from, and an arm that wants an iced start declares a
sub-freezing `sst_pole_k` or a non-zero `sea_ice_thickness_m`. `xclsst` and
`yclsst` are both declared at `-999.` explicitly (`icemod.f90:174`,
`oceanmod.f90:117`) and stay there on purpose: the sentinel is now what
`nfluko` and the flux-correction path test, and `world-4ba` made those refuse
it rather than relax toward it. `make_ice_thickness` survives for the case where
a real code-169 climatology IS supplied, and `iceini` writes a warning naming
the CCM3 relation when it uses it (`icemod.f90:307-310`).

Measured on 2026-08-24, first output bin of `run_2b20e3324bb0`, on the model as
it stood before the repair:

    ts   ocean: min=214.832  max=271.900  mean=258.860
    sic  ocean: min=0.000    max=0.000    mean=0.000
    sit  ocean: NH mean=2.533   SH mean=0.139

The modelled ocean started globally at `TFREEZE`. Every ocean cell carried ice
thickness in an 18:1 north-south ratio, which is Earth's Arctic-Antarctic
contrast imposed on a world with 32 degrees of obliquity and unrelated
geography, while ice COVER was identically zero: `xicecc` initialises to 0 and
`mkicec` can only grow it. Thick ice carrying no albedo is not a physical state,
and neither is the field it replaced.

**What it cost.** The one long run escaped to 294.5 K with 0.07 per cent ice,
so the equilibrium survived. The exposure was `sweep_flux_earth`, which had not
been run: every arm would have cold-started from a fully iced ocean, so the warm
branch could only ever be approached from below and the hysteresis never probed
from the ice-free side, and an arm at 0.90 that failed to melt out would have
read as a snowball bifurcation when it was an initial-condition artifact. Both
fields are written to the restart and read back, so the sentinel and the CCM3
thicknesses persisted through every continuation of a run started that way.
`yclsst` is still emitted as ocean code 990, "clim. sst", and with no
climatology supplied that field is the sentinel and means nothing.

## 2. The derived hyperdiffusion timescales were applied 25 per cent too weak

Mechanism II, in its most expensive instance. `check_consistency.py` computes
`base = math.pi*radius/(ntru*wind)/86400.0` and `run_exoplasim.py` writes that
straight into `TDISSD/Z/T/Q`, in 24-hour days.

`dayseccheck` converts only when the values look like days,
`zmax < day_24hr/mtspd`; at `mtspd = 32` that threshold is 2700 s and the staged
`TDISSZ = 1.1143` is far below it, so the conversion branch is taken and the
value reaches the nondimensionalisation in seconds, correctly. What was wrong
was the nondimensionalisation itself: it divided by `day_24hr` where the model's
unit of time is `sidereal_day/TWOPI`, so every timescale was applied over
`1/rotspd` times its namelist value.

Measured against the T21 rule: tau_vorticity = pi * 7645464 / (21 * 5.94) =
192,552 s, and the model applied 240,678 s, exactly 108000/86400. All four
timescales were 25 per cent too long; the ratios between them survived, so the
SHAPE of the damping was right and only the level was wrong.

**Closed by world-rt1.** `plasim.f90:1985-2002` nondimensionalises `tdissd`,
`tdissz`, `tdisst` and `tdissq` as `sidereal_day/(TWOPI*tdiss)`, and `restim`
and `dampsp` take the same treatment at `:1969` and `:1751-1753`. What
`check_consistency.py` verifies is now what the integration applies, and
`config/planet.yaml`'s claim that the rule is planet-agnostic by construction --
radius and wind, no Earth days -- holds of the model and not only of the
derivation.

The confinement requirement was never in danger either way: solving
`x = (dt/(2*kappa*tau))^(1/(gamma-1)) >= 0.6` needs gamma >= 14.78 at the
nominal tau and gamma >= 15.22 at the tau that was actually applied, and
`filter_power: 16` clears both.

The same defect reached `tfrc` and `taucool`. `tfrc` is live and is now declared;
see finding 20. `taucool` is a physical-seconds conversion at
`plasim.f90:4348`, `taucool*day_24hr`, which is the genuine unit role and is
right; it is reached only at `nstratosponge > 0`.

## 3. The ocean's horizontal diffusion ran on Earth's radius, and the arms that ran are mislabelled

`oceanmod.f90` carried `parameter(PLARAD=6.371E6)`, with the comment "Earth
radius (m)", used in `hdiffo` as `zfac=hdiffk(jlev)/plarad/plarad`. `oceanmod`
did `use resmod` rather than `use pumamod`, so it never saw the `plarad` that
`planet_namelist` sets, and being a `parameter` no namelist could move it.

Live and already exercised: `NHDIFF` and `HDIFFK` are written unconditionally by
`run_exoplasim.py`, and the CLIM-16 and CLIM-19 arms ran at 300, 1000 and 3000.

**Fixed by world-mll.** The parameter is gone. `hdiffo` takes the radius from
`planet_nl` by way of `pumamod` (`oceanmod.f90:1376`, applied at `:1394`),
`oceanini` aborts if `nhdiff > 0` with no radius and logs the radius it will use
(`oceanmod.f90:291-296`), and `hdiffk` now means what it says. The grid there is
still angular, so the radius is what supplies the metric.

**What it cost, and what is still mislabelled.** The divisor was `hdiffo`'s own
`parameter(PLARAD=6.371E6)` and the metric that belonged there is this planet's
`PLARAD = 7645464.0`, so the factor is (7645464/6371000)^2 = 1.4401. Every arm
applied 1.44 times the tendency of its nominal diffusivity, so the CLIM-16
bracket was really 432/1440/4320 and the predicted rms heating table in
`exoplasim/notes/forcing-bundle-predictions.md` understates the model by the
same factor, 0.75/2.51/7.53 becoming 1.08/3.61/10.84. Those labels are wrong
wherever that bracket is quoted; world-vho is the correction, and clim-65's arms
have to be declared against the fixed meaning.

This also resolves an open discrepancy. `predict_ocean_terms.py:99` builds its
operator with `a = radius_earth * EARTH_RADIUS_M`, which is Vesper's radius,
while claiming at line 58 to reproduce `hdiffo`'s discretisation, and its P2
eigenfunction check at line 80 validates the offline code against the analytic
Laplacian on that same radius -- so it structurally cannot catch this. The
prediction was +0.03 to +0.12 K and the measurement +0.13 to +0.19 K; dividing
the measurement by 1.44 lands it on the prediction. The 35 to 37 per cent
inter-realisation spread makes that a consistent explanation rather than a clean
reconciliation, but the factor is not in doubt.

## 4. A 5 m snow cap structurally neutered the glacier module that is switched on

`dsmax = 5.00` m water equivalent, applied in `landstep`. The cap is applied to
`dsnowz`, which is the same accumulator `glaciermod`'s `oroini` builds ice
orography from, so `glacieroro` could not exceed 5.0/0.85, about 5.9 m of ice.
The module's stated purpose, "capturing the effects of large ice sheets on
atmospheric circulation" (`glaciermod.f90:1-8`), was unreachable, and
`where (dsnowz(:) > 30.0) dglac = 1.0` at `glaciermod.f90:266` was dead code.

5 m is PlaSim's Earth seasonal-snow ceiling, chosen for a model where perennial
ice arrives as a prescribed mask rather than being grown by this module.
`NGLACIER = 1` and `GLACELIM = 2.0` in the staged namelist, so the module is on.

**Closed by world-cwc, and the cap is lifted rather than merely settable.**
`dsmax` is a `landmod_nl` key (`landmod.f90:107`), broadcast because it is
threadprivate, and `config/planet.yaml` declares
`surface.glaciers.max_snow_depth_m: -1.0`, which `landmod.f90:980` reads as no
limit. The `dsnowz` accumulator is unbounded, `glacieroro` can grow an ice
sheet, and the 30 m branch in `glaciermod` is live code.

**What it cost.** Everything above the cap was diagnosed as snowmelt and pushed
into `dwater`, so into soil water and runoff. Any cell that would have glaciated
instead manufactured perpetual meltwater at up to the full snowfall rate,
inflating runoff and latent cooling and suppressing the ice-albedo feedback
exactly at the margin. This was the finding that most limited what a run on the
pre-fix model could show, and it bears on `notes/glacier-rough-pass.md`, whose
areas are an upper bound for other reasons as well.

## 5. Cloud liquid water was distributed with an Earth scale height

`zzh(:)=700.*ALOG(1.+dqvi(:))`, the CCM3 liquid-water scale height, a length in
metres fitted on Earth. It is used against mid-layer heights built as
`-dt*gascon/ga*ALOG(...)`, which correctly scale as 1/g. The heights shrink by
0.766 and the scale height did not. Both `700.` and the `0.00021` beside it were
bare literals with no namelist route.

**Closed by world-ofn.** Both are `rainmod_nl` keys: `clwhsc`
(`rainmod.f90:71`) and `clwref` (`:73`), used at `:2031` and `:2034`. `clwhsc`
is not a literal but a DERIVATION: at or below zero it is set in `rainini` to
`700.*(gascon/ga)/(287.0/9.80665)` (`rainmod.f90:187`), so the length carries
this world's scale height and an Earth configuration reproduces 700 exactly.
The default is -1.0, so the derivation is what runs unless a run declares
otherwise, and `rainini` echoes the value it used.

Integrated at 25 kg/m2 precipitable water, the Vesper-over-Earth liquid water
per layer ran 7.2 at sigma 0.038, 4.1 at 0.12, 2.6 at 0.21, 1.9 at 0.32, 1.5 at
0.44, 1.2 at 0.57, 1.0 at 0.70, 0.89 at 0.82, 0.81 at 0.92 and 0.78 at 0.98.
**That table is on the wrong sigma set**, `plasim.f90`'s `neqsig == 0` fallback
rather than the rescaled `NEQSIG = 4` set every run on record uses; see the
correction under finding 20. The two agree from level 3 down and differ by a
third and a fifth at levels 1 and 2, which is where a cloud-water profile is
read, so the top-layer entries are the ones the correction moves.

**What it cost.** Column-integrated liquid water barely moved, 5 to 7 per cent
low, but mid and upper cloud carried 1.5 to 4 times the calibrated water at the
same cloud fraction while the lowest three layers carried 11 to 22 per cent
less. That is the wrong direction on both terms that matter: thicker high cloud
raises the longwave trap, thinner low cloud lowers the shortwave reflection.
`dql` is not diagnostic; it sets shortwave cloud optical depth and longwave
cloud emissivity in `radmod`. The correction applied is exactly the one the
finding named, scaling by `g_earth/g` because liquid water should track vapour
and vapour's geometric scale height is 0.766 here.

## 6. Land longwave surface emissivity is exactly 1.0, and is now declared to be

`lwr` carried `zeps(:) = dls(:) + 0.98*(1.-dls(:))` as a bare literal: the
modelled ocean got 0.98 and the modelled land a perfect blackbody, with no
comment, no unit and no source anywhere, and no namelist key.

**Half settled by world-qvu.** The two numbers are named:
`elwland` and `elwsea` in `radmod_nl` (`radmod.f90:95-96`), applied at `:3212`
as `zeps(:)=elwland*dls(:)+elwsea*(1.-dls(:))`, and `run_exoplasim.py` writes
both unconditionally from `config/planet.yaml`'s `land_longwave_emissivity` and
`sea_longwave_emissivity`. **The declared values are 1.0 and 0.98, the literal's
own**, so no run changed and the physical defect below is untouched: what the
fix bought is that a run records what it emitted with, and that a bracket arm
can move it. The honest answer is a per-cell field derived from the lithology
map, and that is not built.

This is a generic Earth-GCM simplification and it costs more here than on Earth.
`config/planet.yaml` names evaporite and playa clastics as barren classes
covering about a quarter of this land, and salt crust and quartz-rich clastics
run 0.90 to 0.95 broadband in the 8 to 14 micron window.

**What it costs.** At eps = 0.92 and a modelled surface at 300 K the surface
emits 37 W/m2 less and reflects about 26 W/m2 of downwelling that the code
discards outright: `zeps(:)=(1.-zeps(:))*dftd(:,NLEP)` at `radmod.f90:3522` is
identically zero over land while `elwland` is 1.0. Net about 11 W/m2 over that
terrain, which scaled by a quarter of the land fraction is **0.4 to 1.1 W/m2
planetary mean**, bracketed on the emissivity range. That is larger than the
adiabatic term `energy_fixer` was adopted to compensate.

## 7. The snow-over-forest albedo was one spectrum-blind factor applied to both bands

`landini` built the forested-snow endmembers as

    albsmaxf1 = 0.5*albsmax1 ; albsmaxf2 = 0.5*albsmax2
    albsminf1 = 0.75*albsmaxf1 ; albsminf2 = 0.75*albsmaxf2

used in `landstep` to set `dsalb(1,:)` and `dsalb(2,:)`, which are what the
radiation reads at `NSTARTEMP=1`.

The endmembers `albsmax1/2` are derived from the k25v spectrum and are correct.
The canopy masking on top of them was not: 0.5 and 0.75 are Earth-Sun broadband
ratios, from upstream's `albsmaxf=0.4` over `albsmax=0.8`, applied identically
either side of 0.75 micron. This is the class `inherited-earth-constants.md`
finding 1 named for codes 175 and 176, appearing here in the SNOW albedo.

**Closed by world-nfh.** The mixture is written as one, per band:
`albsmaxf_b = forcovmx*albforest(b) + (1-forcovmx)*albsmax_b`, and the same with
`forcovmn` for the aged-snow endmember, with `albforest(2)`, `forcovmx` and
`forcovmn` as `landmod_nl` keys (`landmod.f90:80-82`). The two canopy fractions
are separate keys rather than one because multiple scattering between a bright
snowpack and the canopy masks more of it than the same canopy masks of aged
snow. The defaults solve back to upstream's numbers exactly against an Earth-Sun
canopy endmember of 0.15, so an Earth configuration is unchanged;
`run_exoplasim.py` passes `ALBFOREST` from `config/planet.yaml`'s
`vegetation_albedo_bands`, which `analysis/vegetation_albedo.py` derives from
553 ECOSTRESS leaf spectra re-weighted to this star.

Measured on 2026-08-24 from the staged `orogen_T21_surf_0212.sra`: 2048 cells,
min 0.0000, max 0.4988, nonzero on 997 cells with a mean of 0.3120 over those.
The field is live, not the uniform-land default.

**What it cost.** With `albsmax1 = 0.96819`, `albsmax2 = 0.58529` and band
weights 0.38237/0.61763, anchoring on band 1 with this project's own vegetation
band pair implies a canopy cover of 0.542, which gives a band-2 forested
endmember of 0.390 against the old code's 0.293, so 0.097 too dark at full
canopy. Anchoring on band 2 instead put band 1 about 0.16 too bright; one factor
could not serve both. At the land-mean forest fraction that was about 0.019 too
dark broadband and 0.030 at the maximum, one-signed warm, and the band carrying
the larger half of the error is the one this star puts 0.618 of its surface
shortwave into.

## 8. The critical relative humidity is tuned to an Earth-sized grid box

`rainmod.f90:138`, `rcrit(:)=MAX(0.85,MAX(sigma(:),1.-sigma(:)))`, the threshold
for all non-convective cloud at `rainmod.f90:1997` and `:2003`. It is a subgrid
humidity-variance parameter, so it is a function of grid box AREA, and 0.85 is a
T21-on-Earth tuning. This planet's radius is 1.20 times Earth's, so a T21 box is
1.44 times the area and carries more subgrid variance, which means `rcrit`
should be LOWER.

`RCRIT`, `RCRITMOD` and `RCRITSLOPE` are `rainmod_nl` keys
(`rainmod.f90:74-76`) and this project sets none of them.

**Live, and now labelled where it is set.** world-khn's action on this one was a
declaration and not a change: the block at `rainmod.f90:132-137` says what
`rcrit` encodes, that the by-level modifier normalises by NLEV while the
horizontal assumption has no NLAT term anywhere, and that the value is anchored
to T21. No number moved. The re-derivation for another rung or another radius is
still owed.

**What it costs.** The model produces systematically too little stratiform cloud
at a given grid-mean relative humidity, on the single largest lever this scheme
has over the planet's albedo.

## 9. Precipitation re-evaporation is calibrated to an Earth column and an Earth timestep

`rainmod.f90:47`, `gamma = 0.01`, applied at `rainmod.f90:2255`, `:2266`,
`:2277` and `:2288` as `gamma*(qsat-q)*dsigma(jlev)/deltsec2`. `nevapprec = 1`
by default, so it is on.

Two defects in one constant, and both are still there. The physical content is
the fraction of a layer's sub-saturation that falling precipitation fills in
transit, which is proportional to residence time dz/v_t. Here dz goes as 1/g and
is 0.766 times Earth's, and terminal velocity goes as sqrt(g) at the same air
density is 1.143 times, so residence time is **0.67 times Earth's** and the
Earth tuning over-evaporates by about 1.5 times in every layer. Separately, the
tendency is a per-timestep FRACTION divided by `deltsec2` to make a rate, so the
amount evaporated per step carries no time unit at all and the rate is inversely
proportional to the timestep.

**One thing this finding got wrong, corrected by world-j6v.** `rainini` set
`gamma = 0.007` when `NTRU==42 .and. NLEV==10` and left it at 0.01 otherwise, a
43 per cent step in a precipitation constant taken on a rung change with no
derivation on either side. It cannot have been a truncation compensation,
because what `gamma` is worth depends on the STEP LENGTH and not on the
truncation, and this project gives T21 and T42 the same 45-minute step. The
branch is removed rather than re-keyed, so `gamma` is one number at every rung;
it stays in `rainmod_nl` for a run that wants another value. Every run on record
is T21, where the branch never fired, so 0.01 is unchanged.

**What it costs.** The first half lands on P minus E over land, which is the
numerator of the carve criterion. The second means this term's strength varies
by a factor of two across `resolution_timestep_minutes`, which runs 45.0 at T21,
T42 and T85 and then 30.0 and 22.5 at T127 and T170 -- so a comparison across
the top of the resolution ladder in loop D is comparing two different
re-evaporation strengths.

## 10. Salinity reached the model through four compiled constants, not one

`config/planet.yaml` stated that salinity "reaches the model through one number,
the freezing point of sea water". Four Fortran `parameter`s carried it and no
namelist key existed for three of them.

| constant | what it is | across S = 20 to 40 |
| --- | --- | --- |
| `TFREEZE` | freezing point | handled; -0.04 to -0.36 K |
| `CRHOS = 1030.` | density "at S=34.7", in both `icemod` and `oceanmod` | 1016 to 1032 kg/m3 |
| `CPS = 4180.` | labelled sea water; this was FRESH water's value | 4070 to 3950 |
| `CLFI = 3.28E5` | fusion, depressed by Earth brine content | toward 3.337e5 as it freshens |

**Closed by world-9hb.** All four are `icemod_nl` keys (`icemod.f90:70-72`
beside `TFREEZE`), broadcast, and `oceanmod` no longer carries its own copies:
`oceanini` takes `prhos`, `pcps` and `pclfi` as dummy arguments from `icemod`
(`oceanmod.f90:166-174`, assigned at `:243-244`), so there is one definition and
not two. `run_exoplasim.py` writes `CRHOS`, `CPS` and `CLFI` from the declared
`salinity_psu` and `sea_ice_fusion_j_kg`. `CPS` is now sea water's 3990.34 and
not fresh water's 4180, which is a physics change and moves the mixed-layer heat
capacity in the direction the finding named. `CRHOI = 920.` stays a `parameter`
in both modules: it is pure ice's density and carries no salinity.

**What it cost, in two separable ways.** The mixed-layer heat capacity
`CRHOS*CPS` used at `oceanmod.f90:497` and `:900` was 4.305e6 J/m3/K against sea
water's real 4.10e6 at the declared salinity, so it was 5.1 per cent too large
and the DECLARED `mixed_layer_depth_m: 50.0` was thermally 52.5 m of sea water.
`scripts/error_budget.py` carries that depth as a structural item stated in
metres, and with `CPS` corrected the metres and the thermal meaning agree.

Second, the snow-ice flooding threshold at `icemod.f90:1506-1507` is
`1.E3*xsnow > (CRHOS-CRHOI)*xiced`, a DIFFERENCE of 1030 - 920 = 110, so a 1 per
cent density error is a 10 per cent threshold error. Across the declared bracket
the true margin runs 96 to 112, so -13 to +2 per cent about the compiled value: a
fresher ocean floods and converts snow to ice substantially earlier, trading
insulating snow for conducting ice and accelerating basal growth. That is a
sea-ice feedback CLIM-17's bracket is meant to span, and with `CRHOS` on the
salinity it now can.

## 11. The ozone profile was placed at a fixed geometric height

`radmod.f90:100-101`, `bo3 = 20000.` and `co3 = 5000.`, both in metres, applied
in `mko3` at `radmod.f90:1948-1954`. `mko3` builds its height coordinate
hypsometrically with the model's own `ga`, correctly, and then placed the
profile against a fixed geometric 20 km. This world's scale height is 0.766 of
Earth's, so 20 km is 3.6 scale heights here against 2.7 on Earth.

The modelled photochemical ozone maximum is set by pressure-like conditions --
ultraviolet optical depth and three-body recombination density -- not by
geometric altitude, so holding the height fixed is the wrong invariant. Holding
pressure fixed means `bo3` about 15300 and `co3` about 3830.

**Closed by world-ayx, at the pressure-equivalent values.**
`config/planet.yaml` declares `ozone_height_m: 15311.8` and
`ozone_spread_m: 3828.0`, and `run_exoplasim.py` writes them as `BO3` and `CO3`
into `radmod_namelist` directly, with `expected_namelist_keys` covering both.
That route was chosen over `configure(ozone=dict)` deliberately: the library
path at `__init__.py:2416-2419` writes `BO3`/`CO3` only when `ozone` is a dict,
and `run_exoplasim.py` passes `ozone=bool(atmosphere["ozone"])`, so the boolean
path stays and the two keys are written beside it rather than through it.

**What it cost.** Reproducing the L10 sigma levels and the `mko3` recursion
with a 230 to 288 K profile, the fraction of the modelled ozone column deposited
in the top model layer went from 0.616 at Earth gravity to **0.794** at the
geometric 20 km, and layer 2's share halved from 0.202 to 0.105. That
redistributed shortwave heating rather than changing the column, so the budget
effect is small and consistent with `exoplasim/notes/ozone.md` calling it
second-order for the surface; the top-layer heating rate rose about 29 per cent,
and that is the layer setting the modelled tropopause temperature and static
stability. Distinct from the `a1o3`/`aco3`/`toffo3` shape issue that note
records, which is still open: this one was about gravity.

## 12. The sea-ice albedo ramp kept an Earth broadband slope under star-derived endpoints

`seaini`'s two albedo blocks read

    dsalb(b,:) = doceanalb(b)*(1-dicec) + dicec*AMIN1(dicealbmx(b), dicealbmn(b) + 0.025*(273.-dts))

`dicealbmn` and `dicealbmx` are star-weighted per band by `radmod` and
independently verified by `analysis/ice_albedo.py`. That work stopped at the
endpoints. The slope 0.025 per K is Earth's, tuned to Earth's broadband span
`albice - 0.5 = 0.25`, which is exactly 10 K, and neither it nor the 273.0
anchor was a namelist key. `landmod` does the same job correctly with a
FRACTIONAL ramp, which is the structure this should have had.

**Closed by world-cj4.** The ramp is fractional and shared:
`zsicf(:)=AMAX1(0.0,AMIN1(1.0,(TMELT-dts(:))/dicealbdt))` at `seamod.f90:181`,
applied to both bands and to the non-spectral `albice` at `:193-197` and
`:328-330`. `dicealbdt` is the ramp WIDTH and a `seamod_nl` key
(`seamod.f90:37`), so both bands saturate at the same temperature and the warm
end is anchored on `TMELT` rather than on a bare 273.0. The default width is
10.0 K, upstream's implied width at Earth's broadband span, so an Earth
configuration keeps the same saturation temperature.

**What it cost.** Band 1 spans 0.6203 to 0.8702, a span of 0.250, and saturated
10.0 K below freezing; band 2 spans 0.3395 to 0.4729, a span of 0.133, and
saturated 5.3 K below. The two bands therefore reached maximum albedo 4.7 K
apart, a gap set entirely by how far the re-weighting stretched each endpoint and
by nothing about the modelled ice. At 268 K band 2 read 0.4645 against a
consistently scaled 0.4062, so +0.058 in band 2 and about +0.036 broadband. Near
zero in the global mean, because most modelled ice is colder than either
saturation point, and several W/m2 locally in the 263 to 273 K band where the
ice-albedo feedback lives. One tell that 273.0 was a magic number rather than a
physical anchor: it was neither `TMELT = 273.16` nor `TFREEZE`, so at the melting
point the formula returned 0.004 BELOW its own declared minimum.

## 13. The free-convection ocean transfer coefficient embedded Earth's gravity

The Miller, Beljaars and Palmer (1992) unstable-over-ocean branch read

    zrifh(jhor)=(1.+(0.0016*zdth**(1./3.)/SQRT(zabsu2(jhor))/zkblnz2)**1.25)**0.8

In the free-convection limit this collapses to a transfer velocity of
`0.0016 * dtheta**(1/3)`, the classical four-thirds flux law, in which the
0.0016 IS the `(g/theta)^(1/3)` group evaluated at 9.81. This world needs
`0.0016*(12.81/9.81)^(1/3) = 0.00175`, **9.3 per cent higher**. Unlike the
Charnock case there is no compensating length scale, because the four-thirds law
carries no boundary-layer depth.

**Closed by world-e2k.** `freeconv` is a module variable set once in `fluxini`
from this world's gravity, `freeconv = 0.0016 * (ga / 9.80665)**(1./3.)`
(`fluxmod.f90:137`), echoed to the diag, and used in place of the literal at
`fluxmod.f90:308`. It is derived rather than a key, which is the right shape: it
is not a free parameter but a group whose value follows from `ga`.

**What it cost.** `dtransh` feeds both `mkshfl` and `mkevap`, so this was
modelled ocean sensible heat and evaporation together. At the lowest level
height of about 111 m and a sea roughness of order 3e-5, the neutral coefficient
is about 7e-4, and at 2 m/s with a 1 K air-sea difference the free-convection
term already nearly doubles the transfer coefficient; the sensitivity there is
about 0.54 and rises to 1.0 in the fully calm limit. So a 5 to 9 per cent low
bias in latent and sensible heat over the calm tropics.

Structural, unchanged, and worth recording beside it: there is no
free-convection branch over land at all, so calm daytime fluxes over the playa
and evaporite quarter get no gustiness enhancement.

## 14. One Earth-average soil thermal pair, and no per-cell field behind it

`landmod.f90:143`, `soildiff = 1.8` W/m/K, and `landmod.f90:146`,
`soilcap = 2.4E6` J/m3/K, applied uniformly at `landmod.f90:835-836`. These are
moist mineral soil, Earth's global average, on a project that already derives
albedo and roughness per cell from lithology.

**Half settled, and world-sy9 is OPEN.** The reachability half is done: all six
thermal constants are `landmod_nl` keys and broadcast --
`soildiff`, `sicediff`, `snowdiff`, `soilcap`, `sicecap` and `snowcap` at
`landmod.f90:143-148`, declared at `:126-131`, in the namelist at `:251-253`,
broadcast at `:376-381` -- at the values the literals carried, so no run
changed. `rhosnow` joined them under world-1pl; see finding 25.

**What is still open is the physics, and it is the whole of the finding.** These
are scalars and cannot be fields, so one value covers every lithology whatever a
run declares. A per-cell thermal inertia derived from the lithology map is a
source change to `landmod`, not a namelist value, and nothing has been built.

**What it costs, today.** Thermal inertia I = sqrt(k*rho*c) gives 2078
J/m2/K/s^0.5 for the default against roughly 625 for a dry playa or salt crust,
so about a quarter of the modelled land carries **3.3 times too much thermal
inertia** and its diurnal and seasonal surface temperature swing is damped by
roughly that factor. That feeds evaporation, which feeds P minus E, and it feeds
the dust emission threshold.

---

# Live, bounded or diagnostic

## 15. The archived surface albedo was the Earth-Sun broadband value, not the one the radiation used

`landmod` declares non-spectral `albsmin`, `albsmax`, `albsminf`, `albsmaxf`,
`albgmin` and `albgmax` beside the banded copies. `landini` re-derives only the
banded copies from the spectrum; the six scalars are never touched. `landstep`
computed `dalb(jhor)` from them, and since `radstep` runs before `surfstep` in
the timestep, radmod's correct overwrite of `dalb` was discarded by landmod
every step. The landmod value is what `outmod` accumulates and writes as code
175 (`outmod.f90:662`, `:696`, `:1449`, `:1874`).

Radiation itself reads `dsalb` and was unaffected, so this was archival only.
The stored field read +0.068 over cold snow, +0.057 over snowy forest and +0.061
to +0.068 over cold and melting glacier, one-signed and only over the
cryosphere.

**Closed by world-e5p.** Under `two_band_albedo` `landstep` sets the archived
`dalb` from the bands the radiation actually used,
`dalb(:)=zsolars(1)*dsalb(1,:)+zsolars(2)*dsalb(2,:)` (`landmod.f90:567` and
`:746`), so code 175 carries the flux-weighted combination of the fields that
were integrated with rather than a parallel broadband set.

Code 175 still collides with `dalbcl1` on input (`surfmod.f90:278`), so the same
code means two different things in and out. That is recorded at the collision
now (`surfmod.f90:265-277`) rather than only here: an output 175 is the
all-surfaces albedo including snow, forest, glacier and sea ice, and an input
175 is a background albedo that the model would then put snow back on top of.

## 16. Mean sea-level pressure is reduced with Earth's standard-atmosphere lapse rate, and is no longer written

`vendor/exoplasim/exoplasim/pyburn.py:259`, `RLAPSE = 0.0065`, used at `:1873`
and `:1876` and mirrored at `:2713` and `:2716`, with the ECMWF cold-surface
correction `tstar[tstar<255.0] = 0.5*(255+tstar[tstar<255.0])` beside it.
Everything else in that block, `gascon` and `gravity`, is planet-aware and
passed down from `configure()`; `RLAPSE` is the sole exception in the whole
postprocessor, and `pyburn.postprocess` has no lapse-rate argument. All of that
is still true of the postprocessor.

This world's dry adiabat is 12.75 K/km, so the same fraction of it that 6.5 K/km
is of Earth's gives about 8.5 K/km. The reduction exponent `alpha = gascon*RLAPSE/gravity`
goes from a consistent 0.1904 to the coded 0.1456, which is +5.3 hPa at 2 km of
terrain and +10.3 hPa at 4 km, one-signed and orography-shaped.

**Closed by world-ld1, by not writing the field.** Code 151 is deliberately
absent from both `REGULAR_CODES` and `SNAPSHOT_CODES` in `run_exoplasim.py`,
with the reason recorded at the list. Nothing in the project read `psl`, and a
field that is written, carries a CF standard name and looks authoritative is a
trap whether or not anyone has walked into it yet. `RLAPSE` stays wrong in
`pyburn` and becomes live again the moment anyone adds 151 back or calls
`pyburn.postprocess` for it directly. Same class as
`inherited-earth-constants.md` finding 3, which found 6.5 K/km in three project
scripts and not in the postprocessor; world-ld1 named the two Earth constants
`psl` rests on where they sit.

## 17. Earth's lapse rate and a 288 K surface set every cold start's initial profile

`setzt` builds the restoration and initial temperature profile from `tgr`,
`ALR` and `dtrop`. All three are namelist-reachable, through `plasim_nl` for the
first and third and `planet_nl` for `ALR`, and this project used to set none of
them, so all three arrived at ExoPlaSim's Earth defaults of 288 K, 6.5 K/km and
12 km. `setzt` is reached through `initfd`, which is gated on `nrestart == 0`,
so this is the cold-start path only.

**Closed by world-wmw.** `config/planet.yaml` carries a `cold_start_profile`
block and `run_exoplasim.py` writes `ALR` into `planet_namelist` and `TGR` and
`DTROP` into `plasim_namelist` from it, with all three in
`expected_namelist_keys` so a continuation cannot drop them. `TGR` reverting
matters beyond the profile, which is why it is under the gate: see the second
use below.

`restim = 0` so `damp = 0` and `sr` is an initial state rather than a
restoration; `dtep = dtns = 0` so the latitudinal structure is null. This
world's dry adiabat is 12.76 K/km against Earth's 9.76 and its scale height is
0.77 times, so a 6.5 K/km profile with a 12 km tropopause was a poorer first
guess here than it is on Earth. It washes out, so this cost spin-up rather than
an answer -- but `config/planet.yaml` records that the resolution-ladder
stability boundary was measured entirely on cold starts, and that boundary is
the ceiling for exactly this transient.

Second and smaller use, and the one that is not confined to cold starts:
`plasim.f90` does `zdlnp = zmeanoro/gascon/tgr; psurf = exp(log(psurf)-zdlnp)`,
so the declared 1 bar is reduced for mean orography using `TGR`. At a mean
orographic geopotential of order 2500 m2/s2 that is about 200 Pa, 0.2 per cent
of surface pressure, but it runs at every start, restart included, so the
realised global-mean surface pressure depended on an undeclared Earth number.

## 18. The calendar was a 360-day Earth year that `calini` never updated

`calmod` declares its own `n_days_per_month`, `n_days_per_year`,
`m_days_per_year`, `m_days_per_month` and `mtspd` beside `pumamod`'s, with the
same names. `calini` copied two of the five and derived `mtspd` from the pair it
had not copied, so the calendar's own trio stayed at Earth's 360 days, 30-day
months and an `mtspd` computed against them, while `readnl` had set `pumamod`'s
to this world's.

What that cost, measured before the repair. `mtspd` in calmod was `146*40/360`,
truncating to 16, so the calendar year was `360*16 = 5760` steps against the
orbit's 5850: 1.54 per cent per orbit, a full year out after about 65 orbits.
`tcalday` came out at 43200 s, so a "calendar day" was half a 24-hour day and
`step2cal30` reported hours 0 to 11. `cal2step` and `step2cal30` were not
inverses -- one encoded with `mtspd = 32`, 146 and 15, the other decoded with
16, 360 and 30 -- so year 1 month 1 day 1 encoded to step 4672 and decoded as
`23-Oct-0001`, and `outmod` packed `ihead(3)` from that while `ihead(8)`
reported `pumamod`'s 183, one eight-word header carrying two days-per-year.

`pumamod`'s `m_days_per_month` had no setter anywhere and stayed at 30 against a
183-day year, so twelve months did not span the year in either module.

**The repair.** `calini` copies all five and takes `mtspd` from the caller
rather than deriving it, so the calendar day IS the 24-hour day, `tcalday` comes
out at `day_24hr` by construction, and the calendar year is
`m_days_per_year * mtspd`. `readnl` sets `m_days_per_month` as the year over
twelve ROUNDED UP, so twelve months always cover the orbit and the last is the
short one; rounding down gives a thirteenth month. `cal2step`'s simplified
branch uses `m_days_*` and takes `kyea-1`, which makes it the exact inverse of
`step2cal30` -- verified over three years of steps -- and makes year 1 month 1
day 1 encode to step 0.

That last one was not only a labelling defect. `plasim`'s cold start sets
`nstep = n_start_step` from that call and `radmod` takes the orbital phase from
`mod(nstep,n_steps_per_year)`, so every cold start began at phase 0.799, four
fifths of an orbit past the `meananomaly0` the config declares. It now begins at
phase 0. This changes the initial condition of every new cold start; by rule 7
every existing build and run is disposable and nothing is owed.

A residual remains and is irreducible: the calendar year is 5856 steps against
the orbit's 5850, 0.10 per cent, because `m_days_per_year` and `mtspd` are
integers and the orbit is 182.5 24-hour days. It was 1.54 per cent.

**What reads it.** `pyburn` takes the netCDF time axis from `header[6]`, the raw
step counter, and radiation takes its orbital phase from `n_steps_per_year`
independently, so no published number went through the broken calendar. It
becomes load-bearing the moment anyone enables climatological ozone, t-nudging
or flux correction, the land surface annual cycle, or prescribed ice and SST --
all of which interpolate through `step2cal30`.

The related trap is closed rather than recorded: `n_days_per_year == 365`
switched the whole module to Earth's Gregorian calendar -- `mondays`, the
400/100/4 leap rule, `Jan..Dec` -- and `N_DAYS_PER_YEAR` is derived from this
world's flux and rotation period, so nothing structurally kept it off 365 and
there was no guard and no message. `calini` now aborts on it and says why.

## 19. The glacier persistence test was one run segment, which here is half an Earth year

`glaciermod` initialised `persistflag(NHOR) = .TRUE.`, `glacierstep` only ever
cleared it, and `glacierstop` converted whatever survived to `dglac = 1.0`.
**`persistflag` was not written to the restart**, so the clock reset to a clean
slate at every segment boundary.

The module header says the criterion is "some snow cover continuously for an
entire year" (`glaciermod.f90:10-11`). What it was, was continuously for one
INVOCATION of the model, so the duration was the caller's segment length and was
reachable from no namelist.

**Closed by world-qpe.** `persistt(NHOR)` counts SECONDS of continuous cover at
or above `glacelim` (`glaciermod.f90:64`), incremented by `deltsec` and reset to
zero when the snowpack drops below the threshold or the cell stops being land
(`:656-663`, `:268`); it is a restart record, read at `:163`, so the clock
survives a segment boundary. The duration is the `glacier_nl` key `glacpersist`,
in ORBITS, converted once in `glacierprep` to `glacpersec = glacpersist *
m_days_per_year * day_24hr` (`:120`) and refused if non-positive.
`config/planet.yaml` declares `surface.glaciers.persistence_orbits: 1.0`, and
`run_exoplasim.py` writes `GLACPERSIST` unconditionally so the staged namelist
records the criterion the segment ran with.

The threshold `GLACELIM` is reachable and is set to 2.0, which is ExoPlaSim's
Earth snowline default restated. That half is unchanged: the DEPTH is still an
Earth number, and only the DURATION was fixed.

## 20. An undeclared Rayleigh sponge of 20 and 100 days on the top two layers

`plasim.f90:1538-1540`, run before the namelist read:

    tfrc(1) = 20.0*day_24hr*frcmod
    tfrc(2) = 100.0*day_24hr*frcmod

`TFRC` and `FRCMOD` are in `plasim_nl` and this project set neither. Note the
asymmetry: the NLEV=20 equivalent at `:1561-1570` requires `nrdrag = 1` to
switch on, while the NLEV=10 case, which is the configuration this project runs,
is unconditional.

An inherited Earth-tuned top-of-model sponge, live in every run, on the layers
at sigma 0.025 and 0.094 -- 2.5 and 9.4 kPa on a 1 bar surface -- appearing in
none of this project's declared parameters, unlike the filter, the
hyperdiffusion and the timestep, which are all argued in `config/planet.yaml`.
Whether 20 days suits a 30-hour rotator at 1.31 g is the kind of question this
project's conventions say is answered rather than inherited.

**Closed by world-aee, and the unit was decided rather than inherited.**
`config/planet.yaml` declares `rayleigh_sponge_rotations`, one entry per level
from the top down and zero for no drag, and `run_exoplasim.py` writes the whole
per-level list into `TFRC`. The declared quantity is ROTATIONS and not days,
because the wave field the sponge absorbs has its timescale set by the
rotation; the values are 20.0 and 100.0 rotations, which reproduce upstream's
numbers on a 24-hour rotator and keep the sponge at the same number of rotations
on this one. Keeping them in days would have tightened the sponge to 16
rotations. Mechanism II's error reached `tfrc` too and world-rt1 fixed it there
with everything else.

**THE SIGMA SET IN THIS DOCUMENT WAS THE WRONG ONE.** Measured 2026-08-24
against `run_2b20e3324bb0/MOST_DIAG.00001`, whose vertical table reads 0.02500,
0.09372, 0.18812, 0.29717, 0.42058, 0.55432, 0.69069, 0.81826, 0.92191, 0.98282,
and whose namelist carries `NEQSIG = 4` and `PTOP = 5000.0`. The set this
document used -- 0.038, 0.12, 0.21, 0.32, 0.44, 0.57, 0.70, 0.82, 0.92, 0.98 --
is `plasim.f90`'s `neqsig == 0` fallback, which runs the same polynomial
UNRESCALED and puts the model top at 7660 Pa instead of 5000. Every run on
record is on the rescaled set. The two agree from level 3 down and differ by a
third and a fifth at levels 1 and 2, which is exactly where a top-of-model
sponge and a cloud-water profile are read. Finding 15's per-layer liquid-water
table is computed on the fallback set and is `world-ofn`'s to redo.

## 21. The boundary-layer mixing length is a fixed number of metres in a sigma scheme

`fluxmod.f90:47`, `vdiff_lamm = 160.`, the asymptotic Blackadar mixing length in
metres, used at `:823` and `:846` with `ztscal` at `:61`:

    zzlev = -gascon*ztscal*ALOG(sigmah(jlev))/ga
    zmixm = vdiff_lamm*vonkarman*zzlev/(vdiff_lamm+vonkarman*zzlev)

`zzlev` scales as 1/g and `vdiff_lamm` does not, so the surface-layer to
free-atmosphere crossover sits at a fixed geometric 400 m rather than at a fixed
fraction of a boundary layer that is 0.77 times as deep here. Diffusivity goes
as the square of the length: at the lowest interface, sigmah about 0.966, the
length ratio is 0.845 so K is 71 per cent of the Earth-equivalent, and by
sigmah 0.7 the deficit has fallen to 8 per cent. Under-mixing concentrated
exactly where surface flux is redistributed.

**Live, and named by world-e2k.** `vdiff_lamm` and `ztscal` are both `fluxmod_nl`
keys and broadcast, and the declaration at `fluxmod.f90:40-46` records that 160
m is ECHAM's fit to Earth's free troposphere, that it is not obviously
transferable, and that it is left at ECHAM's value because nothing in this
project has measured a replacement -- not because 160 m has been shown to apply
here. No number moved. The same batch derived `freeconv` from `ga`, which is
finding 13; the distinction is that the four-thirds group FOLLOWS from gravity
while a mixing length does not, so one was derived and the other declared.

## 22. Snow masks the surface at 0.01 m water equivalent regardless of what it buries

`landmod.f90:505-510` and `:672-677`,
`dalb = dalbclim + (zalbsnow - dalbclim)*dsnow/(dsnow+snowcovz)`. Half the snow
albedo effect at 0.01 m water equivalent, about 3 cm of modelled snow at
`rhosnow = 330`. This is a snow-cover-fraction parameterisation whose scale
should track the roughness of what is being buried, and this project's own
staged `orogen_T21_surf_0173.sra` gives land roughness from 0.0038 m to 7.57 m
with a land mean of 2.0 m.

A single 3 cm masking depth turns a 7.6 m rough cell half-white on the first
snowfall. One-signed, too reflective too early, and it applies to the spectral
`dsalb` as well as to `dalb`.

**Live, and named by world-khn.** The literal was at seven sites and is now the
`landmod_nl` key `snowcovz` (`landmod.f90:25`), broadcast, declared at `:17-24`
as the depth at which half the cell is covered and therefore a function of cell
size, anchored to T21 with no NLAT term here or upstream. The value is unchanged
at 0.01 and nothing in this project sets it, so the physical defect above is
exactly as it was; what the fix bought is one definition and an arm.

## 23. Convective cloud cover is an Earth mm/day regression evaluated in 30-hour days

`rainmod.f90:1952-1954` with `zcca=0.245`, `zccb=0.125` at `:1931`:

    zrfac = solar_day * 1000.0 ! convert m/s into mm/day
    zcctot(:)=zcca+zccb*log(dprc(:)*zrfac)

Credit where due: the fork replaced a hardcoded 86400 with `solar_day` here. The
three uncompiled convection variants that still carried `zrfac=86400.*1000.` --
`rainmod_bm.f90`, `rainmod_mca.f90` and `rainmod_kuo_old.f90` -- were deleted
under world-cmz, so the fork has one convection scheme and one day length. But
the regression was fitted on Earth precipitation rates in Earth days, so neither
day length is defensible on the fit's own terms.

Using this world's solar day inflates the log argument by 1.26 and adds a
uniform **+0.029 to convective cloud fraction** wherever the result is not
clamped to [0.05, 0.8]. The point is not which choice is right but that a free
choice here is worth 3 per cent absolute cloud cover.

**Live, and named by world-khn.** `zcca` and `zccb` are still `parameter`s and
no run can move them; the declaration at `rainmod.f90:1926-1930` records the
part that is a resolution problem rather than a day-length one -- the rate the
pair is fitted against is a CELL MEAN, so the same simulated storm spread over a
smaller cell gives a larger rate and more cover, and the cloud fraction of a
convecting region is rung dependent through this pair alone.

## 24. Saturation now branches on ice below the melting point

The model carried one Magnus-Teten triple, `ra1`, `ra2`, `ra4`, the coefficients
over LIQUID water, and used it at every saturation site in `rainmod`, `landmod`,
`seamod`, `simba` and `fluxmod`, while the latent heat already switched to `ALS`
below `TMELT` at four sites in `rainmod` and one in `fluxmod`, with the
Clausius-Clapeyron derivative beside them the liquid one multiplied by L_s/cp.
Saturation over ice is about 25 per cent below saturation over liquid at 250 K,
so cold-cloud condensation was systematically over-produced and the
thermodynamics disagreed with its own latent heat.

`pumamod` now carries `ra1i`, `ra2i`, `ra4i` beside the liquid triple, in
`planet_nl` on `p_earth.f90`, at the standard over-ice 610.66/21.875/7.65.
That triple was in the tree already: `p_mars.f90` carried it in the LIQUID
slots, which is how the set was identified. `p_mars.f90` and `p_exo.f90` have
since been deleted under world-58v and world-cmz -- neither was ever a legal
build value -- so `p_earth.f90` is the only planet module and the only place
the pair of triples has to exist.

The selection goes through three elemental functions in `pumamod`, `ra1s(T)`,
`ra2s(T)` and `ra4s(T)`, so the branch is per gridpoint and every call site
already evaluates an exponential beside it. Thirty-two sites in `rainmod`, two
in `landmod` and one in `simba` take the temperature already in the expression.

Two families deliberately do NOT branch on temperature, and that is what
"consistent with the latent heat already in use" means here:

- `fluxmod`'s two sites take the phase from the `where` arm they are in. That
  mask is `dt > TMELT .or. dls < 0.5`, so a cell with `dls < 0.5` takes the
  liquid arm however cold it is, because the water under a partial ice cover is
  still water and `icemod` owns the ice surface. Branching on temperature there
  would put an ice saturation under a liquid latent heat on exactly those cells.
- `seamod`'s two sea-surface saturations stay liquid, for the same reason.

The audit's earlier reading that this "cannot be fixed from config" was right
about config and beside the point: it needs a second branch rather than a
different single set, so it was always a source change. `alv`, `als` and `tmelt`
are in `p_earth`'s `planet_nl` now, put there by world-58v, so the latent heats
and the switching temperature are settable beside the coefficients.

This is a physics change and it moves the climate: cold-cloud condensation falls
where it was over-produced. A correct term that worsens an agreement is
information.

## 25. Snow and glacier densities are Earth compaction values

`landmod.f90:142`, `rhosnow = 330` kg/m3, and `glaciermod.f90:78`,
`rhoglac = 850` kg/m3. `rhosnow` converts water equivalent to physical snow
thickness for the top-layer heat capacity blend (`landmod.f90:872`, `:879`,
`:922`, `:929`), so it decides how much the modelled snowpack insulates the
soil; `rhoglac` converts water equivalent to ice thickness for the orography
(`glaciermod.f90:422`). Both are settled densities set by overburden compaction,
which scales with gravity, and this world's is 1.306 times Earth's, so settled
snow would be denser and firn would densify faster.

The fork gets the GRAVITY right in `oroini`, with a comment explaining why, and
then multiplies it by an Earth-compacted density. Worth roughly 10 to 20 per
cent on snow insulation thickness and on glacier relief.

**Live, and named by world-1pl.** `rhosnow` is a `landmod_nl` key and `rhoglac`
a `glacier_nl` key, both broadcast, each with a header block saying what it
converts and that it is Earth's compacted value. Neither number moved and
nothing in this project sets either, so the 10 to 20 per cent stands; a
gravity-scaled compaction density is a derivation nobody has done.

## 26. The lead-closing scale and the thickness clamp were Earth Arctic tunings

`zh0 = 0.5` in `mkicec`, credited to Hippler 1979, with the same 0.5 m in
`mkicecf`. It is the whole of the model's lead parameterisation: a cell's
compactness closes with an e-folding of 0.5 m of new ice growth, and `icestep`
then thresholds it at `thicec = 0.5` into a hard mask. So the growth scale alone
sets how much modelled ice must grow before a cell flips to iced for albedo and
roughness, roughly 0.35 m. On a world whose year is half Earth's the ice grown
per season differs, and this Arctic-calibrated growth scale is what converts it
into cover. `thicec` was a namelist key and unset; `zh0` was not a key at all.

Beside it, `xmaxd = 9.0` m was Earth's Arctic multi-year maximum, and its
enforcement was worse than a clamp: the excess was melted and `getiflx` took the
GLOBAL area-weighted sum of that heat and redistributed it onto every other cell
with ice below 9 m, which is a non-local heat transport with no physical
carrier, and one that silently dropped the remainder whenever the demand
exceeded the capacity.

**Closed by world-12c, in both halves.** The lead scale is the `icemod_nl` key
`hlead` (`icemod.f90:113`), declared at `:105-112` and broadcast, and
`run_exoplasim.py` writes it from `config/planet.yaml`'s
`sea_ice_lead_closing_m: 0.5` -- so the value is unchanged and is now declared
Earth Arctic rather than compiled Earth Arctic. The thickness limit pays for
itself LOCALLY: the latent heat comes out of the same cell's conductive flux to
the ocean (`icemod.f90:906-931`), which is conserving and is the only sink the
cell has, and `getiflx`'s global redistribution is gone. And the limit is
switched OFF for this world: `xmaxd` is an `icemod_nl` key and
`config/planet.yaml` declares `sea_ice_max_thickness_m: -1.0`, which
`icemod.f90:920` reads as no limit. That matters beyond the clamp, because a
positive `xmaxd` also zeroes the conductive flux in `mkcflux` and `skintemp` once
ice reaches it, stopping basal growth outright.

Dormant at the climate measured here -- thickness peaked at 4.46 m in spin-up
and 0.78 m at equilibrium, so the branch never fired -- and it was live in any
cold sweep arm, which is what the fix was for.

## 27. Runoff velocity constants were Earth-fitted, with gravity entering through the geopotential

`zcvel = 4.2` and `zcexp = 0.18`, duplicated between `landmod`'s `roffini` and
`glaciermod`'s `oroini`, used as `u = zcvel/zdx * |grad(zoro)|^0.18`. `zoro`
there is GEOPOTENTIAL, not elevation, so the slope term carries `ga`: the same
topographic slope gives 1.306 times the value here, worth 1.306^0.18 = 1.049 on
the velocity. Separately, open-channel velocity scales as sqrt(g), so a
physically consistent `zcvel` is about 1.14 times larger. The pit-filling
increment was likewise 1 m2/s2, so 0.078 m of elevation here against 0.102 m on
Earth.

**Closed by world-529.** The three are `landmod_nl` keys -- `roffvel`,
`roffexp`, `roffpit` (`landmod.f90:98-100`) -- and `glaciermod`'s copy reads
landmod's rather than declaring its own (`glaciermod.f90:382-384`), so there is
one definition and the ordering below no longer decides which numbers win. Each
carries this world's gravity explicitly: `zcvel = roffvel*ga**(0.5-roffexp)`
composes the sqrt(g) of open-channel velocity with the `ga**roffexp` the
geopotential slope already contributes, and `zoroinc = roffpit*ga` puts the
pit-filling increment in metres of elevation. The defaults, `roffvel =
2.022845` per sqrt(g) and `roffpit = 0.101972` m, reproduce upstream's 4.2 and
1 m2/s2 at Earth's gravity exactly.

`oroini` still runs after `roffini` in `surfini`, so `glaciermod`'s assignment
is the one that survives; it now assigns the same numbers landmod did.

## 28. Every run's log identified the planet as Earth

`print_planet` in `p_earth.f90` has no namelist and nothing parses it, so this
was diagnostic only -- but it goes into `plasim_diag` beside numbers that are
correct. Measured on 2026-08-24 from `run_2c42c68fe9ca/MOST_DIAG.00000`:

    YPLANET="Earth"
    *              Simulating: Earth                 *
    *                     Mass  [10^24 kg]    5.9736 *
    *        Equatorial radius        [km] 6378.0000 *
    *              Mean radius        [km] 7645.4640 *
    *          Surface gravity      [m/s2]   12.8100 *
    *              Bond albedo                0.3850 *
    *   Black-body temperature         [K]  247.3000 *

Earth's mass, equatorial and polar radius, density, Bond albedo, blackbody
temperature, perihelion and aphelion were printed directly above this world's
correct mean radius, gravity and irradiance. Only `plarad`, `ga`, `gsol0` and
`sidereal_day` in that table were the real values. Not merely stale: a different
planet, in the artifact a reader consults to confirm what ran.

**Closed by world-1o4.** Every row is now either a namelist value or is DERIVED
from one (`p_earth.f90:128-175`). Mass, volume and density follow from `ga` and
`plarad` through `g = GM/a^2` on a sphere; run on Earth's own `ga` and `plarad`
the derivation returns 5.964, 108.34 and 5505 against the measured 5.972, 108.32
and 5514, and the residual is the oblateness and the rotation this model does
not have, which is the check that it is the right derivation. Dropped rather
than derived, because the model holds nothing they follow from: equatorial and
polar radii and ellipticity, since this planet is a sphere here; and Bond
albedo, black-body temperature, perihelion and aphelion, since the first is an
outcome of the run and the last three need the star's luminosity, which reaches
the model already collapsed into `gsol0`.

`yplanet` still reads "Earth" (`p_earth.f90:21`), and so does the heading the
table prints under. That is the module's identity, which is mechanism I, not a
number.

---

# Dormant, and what each waits on

None of these affects a run today. Each is recorded because it fires on a change
this project already intends to make, and because a dormant Earth constant reads
exactly like a live one to whoever throws the switch.

**Every item below now says so where the switch is thrown.** world-9d1's action
was the same for each: put the trigger and the consequence at the code a person
reads before flipping it, rather than only in this note. What that produced:

- Berger's Milankovitch series, and the `iyrbp = 1950 - n_start_year` remnant
  that feeds it. The remnant carries the note; `radini` reaches the series only
  at `nfixorb == 0` and `run_exoplasim.py` passes `fixedorbit=True` at every
  call.
- `orb_decl`'s vernal equinox, which hardcodes Earth calendar day 80.5 of a
  365-day year through `lambm0 + (calday - ve/365.)*2.*pie`. `calday` is a
  FRACTION of the year in this fork, so `ve/365` is a constant belonging to
  another planet's calendar. Dormant on one keyword: `keplerian=True` routes
  `solang` to `gen_orb_decl`, which takes its phase from `mvelpp` and
  `meananom0r` and is correct and config-reachable.
- `rainmod_bm`'s `ztaud = 7200.` and `ztaus = 14400.`, Earth's canonical
  two-hour deep and four-hour shallow Betts-Miller relaxation, absolute seconds
  with nothing tying them to this world's day, convective depth or timestep.
  `rainmod_bm.f90` was deleted with the other two uncompiled convection variants
  under world-cmz, so this is no longer a switch that can be thrown: swapping
  the convection scheme means bringing a file back from git and re-deriving
  those two timescales for this world, which is the point the finding was
  making.
- The hurricane diagnostics, covered by world-khn's declaration block at
  `hurricanemod`'s thresholds. The positive form is that this diagnostic is not
  meaningful on this world and stays off; if it is ever wanted its thresholds
  are a re-derivation, not a namelist tweak.
- `carbonmod`, dormant at `NCARBON = 0` though `carbonstep` is called every
  timestep. One thing there was a defect and not only a dormancy: `localprecip`
  used `3.154e9`, Earth's seconds per year times 100, two lines from a
  `timeweight` that correctly uses `m_days_per_year`, so the two disagreed by
  the ratio of the years -- about a factor of two here. That is fixed to the
  model's own orbit. `CO2EARTH`, `PEARTH`, `VEARTH`, `RAD_EARTH` beside a
  correct `plarad`, the 288 K reference in `exp(kact*(tsurf - 288.0))` and the
  `tune1`/`tune2` pair commented as fitted to make the global average match are
  all recorded at the module head as what turning `NCARBON` on would need.
- LSG and the coupler, which are not compiled -- `plasim/CMakeLists.txt:183-184`
  links a four-line `src/lsgmod.f90` stub and `cpl_stub.f90` -- and are
  unreachable twice over: `oceanmod.f90:437` aborts unless `n_days_per_year` is
  exactly 360, and that is a count of SIDEREAL DAYS PER ORBIT and is this
  planet's, 146, not a setting. The note sits at the `nlsg > 0` branch
  (`oceanmod.f90:418-434`). The Earth-hardwired bodies it described are gone:
  the whole `lsg/` tree and `plasim/src/cpl.f90` were deleted under world-6ak
  and world-cmz, so what remains is the stub, the refusal, and the record below
  of what turning LSG on would have to be re-derived from.
- `icemod_template.f90`, a stale pre-fork copy in no build that still carried
  `parameter(TFREEZE=271.25)`, the compiled Earth freezing point CLIM-17
  replaced with a namelist key. It was the file a grep for TFREEZE hit first. It
  was deleted under world-cmz, so the grep now lands on `icemod.f90`'s namelist
  key.
- `plasim/bld/make_plasim`, which named `OCEAN=lsgmod` against `lsg/src`, no
  longer exists.

**Dust would have settled about 1790 times too slowly.** `aeromod.f90` declares
`apart = 50e-9` m and `rhop = 1000` kg/m3, a photochemical haze grain at water
density, and both are in `aero_nl`. `enable_dust_emission` wrote every key in
the SOURCE-FIELD provenance, and `build_dust_source_fields.py` builds that dict
with the fourteen `DUST*` keys and neither of these two, so nothing wrote them
and `aero_ini` validates the emission constants and not the grain. Through
`vels`, which is Stokes, the radius enters squared and contributes 1948, the
density 2.602 and the Cunningham factor 0.354 the other way, so 1793 net; and
sedimentation is the only removal term active at `ldepvel = 0` and
`lwetdep = 0`, so the burden would instead have been set by the
timestep-dependent 99-per-cent-per-step bottom-layer scrub. `mmr2n` was off by
the cube.

`enable_dust_emission` now reads the AEROFILE's own sidecar,
`exoplasim/data/dust/vesper_dust_aerosol.provenance.json`, whose
`namelist_values` block carried `APART = 2.2068e-06` and `RHOP = 2600.0`
unread -- the burden-matched radius `dust_aerofile.py` derives for this planet's
gravity -- and writes both into `aero_namelist`, refusing if the sidecar is
missing or either value is not positive. Two provenance files reach that
function and the distinction is load-bearing: the source-field sidecar carries
the emission law, the aerofile sidecar carries the grain.

A new instance of the class in `aerosol-particle-radius.md`, not the audited
one: that audit is about radmod's copy diverging from aeromod's, and this was
about neither copy ever being set.

**SIMBA is Earth-fitted throughout, with two latent traps.** Unreachable at
`NVEG = 0` and one namelist key away. `rlue = 3.4E-10` kg C/J is applied to broadband downward shortwave,
`dfd(jhor,NLEP)` at `simba.f90:450`, rather than to the photosynthetically
active fraction, which under this star is roughly 0.30 to 0.33 against about
0.45 for the Sun, so light-limited modelled GPP would be overstated by about 1.4
times; the module header says so itself at `simba.f90:26-30`. `ct_crit = 5.0` gives a linear ramp from 273.16 K and then
flat forever, with no temperature optimum and no high-temperature cutoff;
`q10 = 2.0` is referenced to 283.16 K; `co2_ref = 360` ppm and `co2_sens = 0.3`
are Harvey (1989)'s Earth beta factor; `zlaimax`, `cveg_k`, `cveg_l`, `cveg_a`
and `cveg_f` are Earth biome fits; `valb_min`, `valb_max` and `vsalb_min` are
Earth-Sun broadband albedos, and `vsalb_min` is a THIRD independent forest-snow
albedo constant inconsistent with both finding 7 and anything derived from this
star. `tau_veg = 10` and `tau_soil = 42` years (`landmod.f90:117-118`) are scaled
correctly and planet-aware at `landmod.f90:396-397`, which is the trap in the
other direction:
it converts Earth-calibrated turnover times into VESPER years, so 10 and 42
become 5 and 21 Earth years of physical turnover. The two latent traps:
`vegstep` (`simba.f90:549`) overwrites `dwmax` every timestep from SIMBA's own
Earth range, silently discarding the pedology-derived code-229 field the whole
bootstrap-to-baseline sequence exists to produce; and `vegstep` sets `dalb`
but not `dsalb`, so under `two_band_albedo: true` the radiation reads `dsalb`
only and SIMBA's vegetation albedo would have no radiative effect at all while
its roughness and bucket would.

**Both traps are disarmed by world-9hv, and the Earth fits are unchanged.**
`vegini` refuses `nveg == 2 .and. nstartemp == 1` outright
(`simba.f90:198-206`), so the incompatible pair cannot run silently; and
`nvegwmax` (`simba.f90:101`) decides who owns `dwmax`, defaulting to SIMBA,
which is upstream's behaviour, with both arms announced at initialisation
(`:213-220`). Setting `NVEGWMAX = 0` keeps the pedology-derived code-229 field.
Everything else in this paragraph -- `rlue` on broadband rather than
photosynthetically active flux, `ct_crit`, `q10`, `co2_ref`, `co2_sens`, the
biome fits and the three Earth-Sun albedo constants -- is unchanged and fires on
`NVEG`.

**Flux correction relaxes toward the constructed Earth ice field.** `nfluko` is
the model's actual q-flux mechanism. `nfluko=1` reads codes 709 and 903, which
do not exist here so they stay zero. `nfluko=2` (`mkflukoi`, `icemod.f90:1735`)
relaxes ice toward `xcliced2`, which is finding 1's CCM3 field, and `addfci`
(`:1758`) compares SST against the interpolation of the `-999` `xclsst`.
`config/planet.yaml` explicitly wants a q-flux bracket and
`analysis/error_budget.json` carries no-q-flux as a structural term, so this is
directly on the path of work already planned.

**The resolution-tuned radiation constants are pinned at T21's row by two
independent mechanisms.** `radmod.f90:968-1030` sets `tswr1`, `tswr2`, `tswr3`
and `th2oc` per (NTRU, NLEV), and every branch is guarded by
`if(NDCYCLE==1) jtune=0`. `ndcycle` defaults to 1 (`radmod.f90:205`) and nothing
sets it, so every branch of that table is dead and the module defaults at
`radmod.f90:72-75` stand at every rung -- and those defaults ARE the T21/L10 row.
Second mechanism: `run_exoplasim.py` writes `TSWR3` as a base times
`cloud_absorption_scale`, hardcoding T21's base. **world-ys9 is OPEN**, and this
is one of the two findings in this document that is not settled. Correct today at T21. At T42 the model would use `tswr1`
0.077 against the table's 0.089, `tswr3` 0.0055 against 0.0048, and `th2oc`
0.024 against 0.0285, so 16 per cent on the longwave water-vapour continuum.
`th2oc` is a new member of the cloud-tuning class of
`inherited-earth-constants.md` finding 2 and is not on the corrected list.

**Berger's Milankovitch series and an Earth day-80.5 equinox.**
`orb_params` (`radmod.f90:3709-4202`) computes Earth's orbital elements from a
year AD, reached from `radini` when `nfixorb == 0`; `run_exoplasim.py` sets
`fixedorbit=True`, so it is unreachable. `orb_decl` (`:4209`) carries
`lambm = lambm0 + (calday - ve/365.)*2.*pie` with `ve = 80.5`, Earth's
Jan-1-to-equinox phase, which does not scale with the model calendar; dormant
because `keplerian=True` routes `solang` at `radmod.f90:1837` to
`gen_orb_decl` (`:3549`), which takes its phase from `mvelpp` and `meananom0r`
and is correct and config-reachable. world-9d1 put the trigger and the
consequence at `radmod.f90:4270-4275`, where the constant is. Both become live
on one keyword.

**Betts-Miller's Earth adjustment timescales.** `rainmod_bm.f90` carried
`ztaud=7200.` and `ztaus=14400.`, the canonical two-hour deep and four-hour
shallow Earth relaxation. The file was deleted under world-cmz along with the
other two uncompiled convection variants, so swapping the convection scheme
means restoring it from git and re-deriving both timescales for this world.

**Hurricane diagnostics are Earth-empirical end to end.** Dormant at
`nstormdiag = 0`, and unchanged. `hurricanemod.f90` declares a wall of Earth
tropical-cyclone thresholds. Two are specifically broken off Earth rather than
merely unfitted: `LAVTHRESH=1.2e-5` per s absolute vorticity (`:64`), against a
planetary vorticity 0.8 times Earth's at 30 hours, so 25 per cent too strict;
and `SIZETHRESH = 30` cells against a T21 cell 1.44 times Earth's area, so a 44
per cent larger storm. The GPI formula carries its Earth normalisations inline.
`RD` is correctly taken from `gascon`, while `CPD=1005.7` at `:51` is hardcoded
and not tied to `acpd`. world-khn's declaration block at `:70-82` records the
four thresholds that are in GRID units rather than physical ones, and world-9d1
put the trigger beside them. The positive form: this diagnostic is not
meaningful on this world and stays off; if it is ever wanted, the thresholds are
a re-derivation and not a namelist tweak.

**LSG and the coupler were Earth-hardwired, were not compiled, and are now
deleted.** `plasim/CMakeLists.txt:183-184` links `src/lsgmod.f90`, a four-line
stub, and `cpl_stub.f90`. The real bodies were `lsg/src/lsgmod.f90`, which
carried `g=9.80665`, `erdrad=6371000.00`, `erdrot=four*pi/86164.`,
`rhonul=1030.`, `tfreez=-1.91`, `cp=4180.` and `entmel=80.*cp`, with a density
polynomial fitted to Earth's T and S range at `sref=35., tref=2.` and a depth
coordinate using `du(k)*rhonul/g`; and `plasim/src/cpl.f90`, which was T21-only
AND Earth-radius-only via `parameter(nxa=64,nya=32)` and
`parameter(radea=6.371E6)`. Both were deleted under world-6ak and world-cmz.
`oceanmod.f90:437` still aborts unless `n_days_per_year == 360`, which this
world's calendar can never satisfy, so the `nlsg > 0` branch remains
unreachable; the list above is what an adopted dynamic ocean would have to be
re-derived from, and `vendor/cgenie` is the candidate this project actually
carries for that role.

**`carbonmod` hardcoded Earth's seconds per year, and still hardcodes Earth's
land fraction.** Unreachable at `NCARBON = 0`, though `carbonstep` is called
every timestep regardless. The seconds-per-year defect was a defect and not only
a dormancy, and world-9d1 fixed it: `carbonmod.f90:198` now reads
`localprecip(:) = (dprc(:) + dprl(:))*m_days_per_year*day_24hr*100.0`, the
model's own orbit, where it used to read `*3.154e9`, Earth's seconds per year
times 100, two lines from a `timeweight` that already used `m_days_per_year`.
The rest is recorded at the module head as what turning `NCARBON` on would need:
`exp(kact*(tsurf - 288.0))` at `:219` takes Earth's global mean as the
activation reference; `:237` divides by `RAD_EARTHSQ*0.29`, where the radius is
correctly the model's but 0.29 is Earth's land fraction against this planet's 50
per cent; `RAD_EARTH=6371220.0` at `:47` sits beside a correct `plarad`; and
`tune1 = 5.41` and `tune2 = 2.20` are commented as "tuning adjustment to make
global average match", so fitted to Earth by construction.

**`icemod_template.f90` was a stale pre-fork copy, and is deleted.** It carried
`parameter(TFREEZE=271.25)`, the `parameter` form from before the fork made it a
namelist key, with no `nseaice`, no threadprivate directives, and no entry in
`CMakeLists.txt`. It was the file a grep for `TFREEZE` hit first and it gave the
pre-CLIM-17 answer. world-cmz deleted it along with the other twelve sources no
configuration compiles; `plasim/bld/make_plasim`, the old unused makefile that
named `OCEAN=lsgmod` against `lsg/src`, went with the build-system removal.

**`newsnow.f90` and `buildice.f90` are orphaned, and no longer hardcode a T21
grid.** Nothing in `vendor/exoplasim` or `exoplasim/` references either; they are
never compiled and never run. Both used to declare `NLAT = 32, NLON = 64`, which
happens to match the current rung and would silently mis-read any other; under
world-zsa both take the grid from `resmod` like every other unit in the tree
(`buildice.f90:4`, `:21-23`, and the same in `newsnow.f90`), so a mismatch is
now a build error rather than a garbled read. `buildice.f90` still adds 400 m
water equivalent to every land cell unconditionally. They are the other side of
the external mass-balance handshake `glacierstop` writes at
`glaciermod.f90:702-709`, which the model produces on every run and nothing
reads. Whether they should exist at all is `clim-53`, which is open; world-cmz
kept them deliberately rather than deleting them with the rest.

---

## Checked and clean

Recorded so they are not re-derived. Several are places where the fork has
already done this work correctly, and knowing which is what makes the findings
above legible as exceptions rather than as the rule.

**Every column amount divides by the model's own gravity.** Ozone
(`radmod.f90:2485`), water vapour (`:2491`), shortwave CO2 (`:2505`), Rayleigh
(`:2514-2515`, `:2736`), longwave H2O/CO2/O3 (`:3226`), CH4/N2O air mass
(`:3238`), cloud liquid water path (`:2426`), and the dust layer thickness.
Rayleigh additionally carries an explicit `9.80665/ga`. Nothing computes a column with
Earth's gravity, and I found no `1/g` omitted anywhere in `rainmod` or
`fluxmod`.

**No Earth boundary field survives into a run.** Codes 169, 210, 211, 709 and
903 all print `Init ... internally` in the diag. `configure()`'s `rm *.sra` does
its job, and the run directory holds exactly the seven project-generated files.
Finding 1 is what the model did INSTEAD, not a surviving file.

**Every derived scaling is planet-general.** `ww = TWOPI/sidereal_day`,
`acpd = gascon/akap`, `adv`, `cv = plarad*ww`, `ct = cv^2/gascon`,
`rdbrv = gascon/RV`, `pnu21` (`plasim.f90:1634-1641`), and `plavor = EZ`, which
is the non-dimensional planetary vorticity in units of `ww` and correct for any
prograde rotator. `akap = 0.286` with `gascon = 287.017` gives
`acpd = 1003.6 J/kg/K`, which is the right answer for this N2/O2/Ar mixture
rather than an inherited error.

**The planetary constants reach the model correctly.** `PLARAD = 7645464.0`,
`GA = 12.81`, `GASCON = 287.017`, `ECCEN = 0.02`, `OBLIQ = 32.0`,
`GSOL0 = 1286.145`, `ROTSPD = 0.8`, `SIDEREAL_YEAR = 15794043.12`,
`NFIXORB = 1`, `NGENKEPLERIAN = 1`, `PSURF = 100000`, all verified in
`run_2c42c68fe9ca`'s staged namelists rather than in documentation. `meananom0`
presets Earth's Jan-1 value at `p_earth.f90:42` and all three drivers pass
`meananomaly0=0.0`, so it never fires.

**Charnock is gravity-correct.** `seamod.f90:320-321` computes
`z0 = charnock*|tau|*gascon*T/(ga*dp)`, exactly `alpha*u_*^2/g` with `ga`
explicit, and `charnock = 0.018` is dimensionless with no g in it. The floor
`dz0sea = 1.5e-5` is a viscous smooth-flow scale, not a gravity quantity.

**The postprocessor is planet-aware everywhere but one constant.** `radius`,
`gravity` and `gascon` are passed explicitly on both branches of
`__init__.py`. Vertical velocity, geopotential height, streamfunction and every
vector transform in `pyburn.py` use them. Finding 16's `RLAPSE` is the sole
exception.

**No Earth calendar reaches any artifact.** The netCDF time axis is
`["time","timestep_of_year","timesteps"]`, raw model timesteps, with no
`calendar` attribute, no epoch and no month names. Seasons are defined by solar
longitude in `analyze_climatology.py` with unequal bin weights handled by
`lib/climatology.py`. Finding 18's calendar defect was confined to the raw
`.srv` headers and the model log.

**The spectral machinery is geometry-free.** `legmod` is fully non-dimensional,
`shtnsmod` carries normalisation constants only and shares `legmod`'s `fsp`
rather than rebuilding it, and `gaussmod` is Newton on the Legendre recurrence.
The physics filter is scale-free in `n/NTRU`, so it damps the same FRACTION of
every rung's spectrum. `tracermod`'s CFL is built from `plarad` and `deltsec`
and is correct for any radius.

**The seasonal cycle is correct and bypasses the broken calendar.**
`zcday = mod(nstep, n_steps_per_year)/real(n_steps_per_year)` with
`n_steps_per_year = 5850`, matching `runsteps_per_orbit` exactly.

**The energy fixer's averaging window is `ntspd` steps**, one solar day of this
planet rather than an Earth day.

**The slab ocean and the sea ice step on the atmosphere's interval, and the
apparent mismatch is a name collision.** Worth recording because the sweep
raised it as a finding and it took source to refute. `icemod.f90:438` reads
`xdt = solar_day / real(ntspd)` and `oceanmod.f90:313` the same for `dtmix`,
which against `pumamod`'s `ntspd = 40` and `solar_day = 108744.83` would give
2718.62 s against the atmosphere's `deltsec = 2700`. It does not, because both
names are rebound on entry: `iceini`'s dummy arguments are `ktspd` and `psolday`
(`icemod.f90:341`), it assigns `ntspd = ktspd` and `solar_day = psolday` at
`:375-376`, and `seamod.f90:133` passes `mtspd` and `day_24hr`. So
`xdt = 86400/32 = 2700 s`, identical to `deltsec = day_24hr/mtspd` at
`plasim.f90:689`. Two different
`ntspd` and two different `solar_day` exist in one model and disagree by
1.0069; reading either pair across the call boundary gives a defect that is not
there.

**The dust emission threshold is already gravity-scaled** by `(g/g_earth)^(1/4)`
from the Shao-Lu threshold minimum, `dustsrc` takes `ga` and `gascon` and builds
the bottom-level height hypsometrically rather than assuming 10 m, and
`rho_a0 = 1.225` is Kok's standardisation density that cancels by construction.
Stokes remains valid at the intended grain size, Re about 6e-4, so the missing
drag-coefficient table costs nothing below about 40 micron.

**Soil layer depths are adequate, checked rather than assumed.** With
kappa = 7.5e-7 m2/s, the annual damping depth is 2.75 m on Earth and 1.94 m
here, so the 12.4 m column is 6.4 Vesper damping depths against 4.5 Earth ones
and the zero-flux bottom boundary is if anything better here. `dztop = 0.20` m
against a 30-hour diurnal damping depth of 0.161 m is 1.25 against Earth's 1.39,
so 11 per cent under-thick, worth a few per cent on the diurnal amplitude.
Reachable via `soildepth` if it ever needs to move.

**The two-band snow, sea-ice, glacier and ground endmembers are derived from the
k25v spectrum at runtime** (`radmod.f90:520-760`) and reach `landmod` through
the namelist, confirmed in the diag echo. The 263.16 K snow aging ramp is
`tmelt - 10` and carries no planet dependence. Findings 7 and 12 concerned what
was applied ON TOP of correct endmembers.

**The compiled convection scheme has no adjustment timescale to be wrong.** Kuo
is a moisture-convergence closure: the moisture increment is accumulated over
`deltsec2` in `mkdca` and redistributed over the same interval in `mkrain`.
There is no CAPE relaxation time, no entrainment rate and no precipitation
efficiency to inherit. There is also no precipitation fall speed anywhere --
`mkrain` moves precipitation through the whole column in one step -- so no Earth
terminal velocity is misapplied.

**The 316 nm shortwave cutoff is reverse-engineered from the Sun and is
nonetheless harmless.** `radmod.f90:218`, `minwavel = 316.036116751`, applied at `:469` and matching
the `zsolar1=0.517` at Teff=5772 K recorded at `:39`, zeroes 0.353 per cent of
k25v's flux. Every absorptance is divided by its band's flux share and
multiplied back, so the partition cancels identically; the residue reaches
`zsolar1` itself at 0.57 per cent relative, worth about 0.06 W/m2 and 0.0002 in
surface albedo.

**The Earth air-mass magnification is harmless.** `zm = 35./SQRT(1.+1224.*mu^2)`
(`radmod.f90:2465`) has `sqrt(1224) = 34.99`, so it reduces to 1/mu everywhere
except below mu about 0.03, where the cell receives under 3 per cent of
normal-incidence flux. Under 0.1 W/m2.

**`tfrc`'s sponge sits on sigma 0.025 and 0.094**, so it is stratospheric
rather than boundary-layer drag, which is what makes finding 20 a declaration
question rather than a boundary-layer one. The claim this entry used to make,
that `tfrc`'s `day_24hr` cancelled correctly, was wrong: it went through the
same nondimensionalisation as the hyperdiffusion timescales and carried the same
`1/rotspd` error, which world-rt1 fixed. `miscmod.f90:75` nondimensionalises
`tnudget` as `TWOPI*tnudget/ww`, which is the same treatment expressed the other
way round and is right.

**Water properties are not Earth properties.** The Magnus coefficients, `RV`,
`ACPV`, `TMELT`, `sicecap`, `sicediff`, `ALS`, `ALV` and
`L_TIMES_RHOH2O` are properties of the condensable. Only their USE is at issue,
in finding 24.

**The shallow-convection sigma diffusivity converts correctly.**
`rainmod.f90:1565` carries the `ga/gascon` factors that make the right Jacobian
for m2/s into sigma2/s, so only the magnitude of `rkshallow = 10.`
(`rainmod.f90:30`, a `rainmod_nl` key) is unjustified rather than the
conversion.

**`t0 = 250.0 K`** (`plasimmod.f90:879`) is the semi-implicit reference
temperature, a reference rather than a constraint, and the `nconvtime`
gravity-wave check derives its limit from it correctly with `gascon` and `akap`.

**`utilities.f90` and `restartmod.f90` contain no planetary constants.**
`makeareas` builds fractional spherical cell areas from latitudes only.

## Two things that are neither findings nor clean

**The dynamical non-dimensional timestep was exact by coincidence, and is now
exact by construction.** `plasim.f90` set `delt = TWOPI/ntspd`, where the
correct value is `deltsec*ww`. `ntspd` counts timesteps per SOLAR day and `ww`
is built from the SIDEREAL one, so the two differ by
`n_days_per_year/(n_days_per_year-1)`, 0.69 per cent here. It did not show,
because `ntspd = nint(solar_day)/nint(mpstep*60)` is an INTEGER division that
truncates to exactly `sidereal_day/deltsec` at every timestep on the rung table
and annuls the error. The truncation stops annulling it above
`n_days_per_year-1` steps per day, which is below dt = 12.4 min, inside the
range the T127 and T170 rungs are headed for. **world-r8o wrote it as
`delt = deltsec*ww`** at `plasim.f90:652` and `:691`, so the identity holds at
every timestep rather than at the ones where an integer division happens to
agree, which is also what `plasimmod.f90`'s declaration of `delt` says it is.
The same lattice is what keeps `solang`'s hour angle advancing exactly 2*pi per
model day with no discontinuity.

**`N_DAYS_PER_YEAR` survives by ordering, and still does.**
`__init__.py:2145-2146` computes `max(int(360.0/rotationperiod/12+0.5),1)*12`,
which is 288 for this world: Earth's 360-day year expressed in Vesper days, not
Vesper's year. It is harmless only because `configure_otherargs` in
`run_exoplasim.py` writes 146 through `otherargs`, which `configure()` applies
last. All three drivers pass that dict, so nobody hits it today.
`NSTPW = int(7200//timestep)` at `__init__.py:2133` is likewise 5 Earth days,
giving 160 steps, which is 4 Vesper days and does not divide 5850;
`run_stellar_cycle.py` works around it explicitly and the other two drivers do
not. The downstream symptom is documented and patched in
`exoplasim/notes/first-output-bin.md`, so this records the Earth-day origin
rather than a new consequence.

## What this audit did not cover

Stated so the coverage is not overread. The `puma`, `cat` and `sam` trees were
not audited, nor the GUI path, the FFT modules, the MPI decomposition, or the
`postprocessor` and `pRT` subtrees beyond the constants they expose. Most of
that is moot rather than outstanding: `puma`, `cat`, `sam`, the GUI sources and
`pRT.py` were deleted under world-6ak and world-cmz as trees nothing builds,
reads or ships, and the MPI build path and the parmode axis went under
world-38b, leaving `mpimod_omp.f90` as the only decomposition. What is still
uncovered is the postprocessor beyond `RLAPSE`. Within the
compiled set, `rainmod`'s cloud-cover diagnostic feeding `dcc` and the leapfrog
`deltsec` against `deltsec2` question were looked at and left; the second is a
water-budget question for `water-and-energy-closure.md` rather than an
Earth-centrism one. Physics that is simply ABSENT was out of scope and is the
subject of `absent-and-inherited-physics.md`.

This audit read source and staged artifacts. It ran no model. Every magnitude
quoted is an estimate from the code and the configuration rather than a measured
perturbation, and the ones most worth measuring are finding 1's sweep behaviour,
finding 2's factor of 1.25, finding 3's factor of 1.44 against the arms already
run, and finding 6's emissivity bracket. None of those measurements has been
made: the fixes below were argued from the code, and every run that would have
supplied a before-and-after is disposable under rule 7 in any case.

## Tasks

Tracked in the `bd` issue tracker, not restated here. Rows carrying the label
`audit:model-earth-centrism` were opened by this audit; the rest were already
open when it ran, from the sweeps of the same day, and this audit corroborates
them rather than duplicating them. Where two rows cover one finding, both are
named because they were filed independently.

The `bd` row is the authority on status; the column below is the state at the
close of the fan-out that worked this audit, and it is here because a reader
arriving at a finding needs to know whether the code still looks like the
finding does. Three rows are open: `world-sy9`, `world-ys9` and `clim-53`.

| finding | id | state |
| --- | --- | --- |
| 0. continuations revert hyperdiffusion and filter power | `world-8bs` | closed |
| 1. the cold start builds a CCM3 Earth Arctic ice field | `world-6fh` | closed; the cold start is declared and ice-free |
| 2. the damping timescales are applied 25 per cent too weak | `world-rt1` | closed |
| 3. the ocean's diffusion runs on Earth's radius | `world-st4`, `world-mll` | closed; arm labels corrected by world-vho |
| 4. a 5 m snow cap neuters the glacier module | `world-cwc` | closed; the cap is lifted, not merely settable |
| 5. cloud liquid water uses an Earth scale height | `world-ofn` | closed; the length is derived from gascon/ga |
| 6. land longwave emissivity is exactly 1.0 | `world-qvu` | named only; the value is still 1.0 |
| 7. the snow-over-forest albedo is one factor for both bands | `world-nfh` | closed |
| 8. the critical relative humidity is an Earth grid-box tuning | `world-khn` | declared only; unchanged |
| 9. precipitation re-evaporation is Earth-calibrated and timestep-dependent | `world-j6v`, `world-khn` | rung branch removed; both physical defects stand |
| 10. salinity is four constants and `CPS` is fresh water's | `world-9hb` | closed; CPS moved to sea water's |
| 11. the ozone profile sits at a fixed geometric height | `world-ayx` | closed; placed on pressure |
| 12. the sea-ice albedo ramp keeps an Earth slope | `world-cj4` | closed |
| 13. the free-convection ocean coefficient embeds Earth's gravity | `world-e2k` | closed; derived from ga |
| 14. one Earth-average soil thermal pair | `world-sy9` | OPEN; settable, still one value for every lithology |
| 15. the archived albedo 175 is the Earth-Sun broadband value | `world-e5p` | closed; the 175/176 collision is documented, not renamed |
| 16. mean sea-level pressure uses Earth's lapse rate | `world-ld1` | closed by not writing 151; RLAPSE unchanged |
| 17. Earth's lapse rate sets the cold-start profile | `world-wmw` | closed |
| 18. the model's calendar is a 360-day Earth year | `world-x1k` | closed |
| 19. glacier persistence is one run segment | `world-qpe` | closed; GLACELIM is still an Earth depth |
| 20. an undeclared Rayleigh sponge of 20 and 100 days | `world-aee` | closed; declared in rotations |
| 21. the mixing length is fixed metres in a sigma scheme | `world-e2k` | declared only; unchanged |
| 22. snow masks the surface at 0.01 m | `world-khn` | declared only; unchanged |
| 23. convective cloud cover is an Earth mm/day regression | `world-khn` | declared only; unchanged |
| 24. no saturation-over-ice branch | `world-ako` | closed; the branch exists |
| 25. snow and glacier densities are Earth compaction values | `world-1pl` | declared only; unchanged |
| 26. the lead-closing scale and the thickness clamp | `world-12c` | closed; the clamp is off and local |
| 27. runoff velocity constants | `world-529` | closed; one definition, gravity carried |
| 28. the run log identifies the planet as Earth | `world-1o4` | closed; every row derived |
| mechanism I, the build compiles `p_earth` | `world-cwu`, `world-58v` | closed; one planet module |
| dust `APART` and `RHOP` were never written | `world-906`, fixed | closed |
| dormant, SIMBA and its two traps | `world-9hv` | closed; both traps refuse or announce |
| dormant, `nfluko` relaxes toward the constructed field | `world-4ba` | closed; nfluko refuses a sentinel |
| dormant, the radiation tuning is pinned at T21 | `world-ys9` | OPEN |
| dormant, the rest, with what each waits on | `world-9d1` | closed |
| dormant, deleting what no build compiles | `world-cmz`, `world-6ak` | closed |
| dormant, `newsnow.f90` and `buildice.f90` | `clim-53`, `world-zsa` | world-zsa closed; clim-53 OPEN |

`world-cwu` recorded `alr` as inert. Finding 17 corrected that, and world-wmw
put `ALR` into the staged `planet_namelist` from the declared cold-start
profile.
