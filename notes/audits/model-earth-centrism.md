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

**Fourteen findings are live and load-bearing, fourteen are live but bounded or
diagnostic, and eleven are dormant traps waiting on a switch this project
already intends to throw.** The most consequential is finding 1, because every
cold start on this world begins from a sea-ice field constructed out of Earth's
Arctic. Finding 0 is not Earth centrism at all and outranks all of them.

The individual constants matter less than the four mechanisms that keep
producing them, and those are stated first because they are what a reader needs
in order to predict where the next one will be.

---

## 0. A converged climatology was integrated on the wrong diffusion operator

Not an Earth-centrism finding. Recorded here because this sweep is what found
it and because it invalidates an artifact the tree currently trusts.

`continue_exoplasim.py` imports and calls five of the six `declare_*` helpers
after `configure()` and omits `declare_hyperdiffusion`. `configure()` rewrites
`plasim_namelist` wholesale on every continuation, so each continued segment
drops back to PlaSim's compiled Earth tuning. The same two drivers also omit
`filterpower`, so `NFILTEREXP` reverts from the derived 16 to the library
default 8.

Checked across all nine run directories on 2026-08-24. One is missing every
key:

    run_2b20e3324bb0   NHDIFF=-    NDEL=-   TDISSZ=-        NFILTEREXP=8
    run_2c42c68fe9ca   NHDIFF=16   NDEL=4   TDISSZ=1.1143   NFILTEREXP=16
    run_52349d4bfb75   NHDIFF=8    NDEL=4   TDISSZ=2.2285   NFILTEREXP=16
    ... six more, all complete ...

`exoplasim/runs/INDEX.json` records that run as `climatology_complete` and
`converged: true`. At T21 the compiled fallback is `nhdiff=15`, `ndel=2`,
`tdissd=0.20`, `tdissz=1.10`, `tdisst=5.60`, `tdissq=0.1`, because the
`NTRU==42` branch at `plasim.f90:1443-1450` does not fire at T21. So that run
integrated with a fourth-order operator instead of an eighth-order one, a
damping cutoff at 0.71 of the truncation instead of 0.38, and humidity damping
7.4 times too strong.

`expected_namelist_keys()` in `run_exoplasim.py:1460-1496` exists for exactly
this class and covers the shortwave gas weights, `O3SCALE`, `TFREEZE`,
`NENERGY`, `NENER3D`, `NCONVTIME`, `NDEALIAS`, `NENERGYFIX` and `PNU`. It does
not cover `NHDIFF`, `NDEL`, `TDISS*`, `NFILTEREXP` or `FILTERKAPPA`, so
`verify_staged_namelists` ran on every continuation and printed "namelists
verified" over both reversions.

## The four mechanisms

**I. The build compiles Earth's planet module.** `plasim/CMakeLists.txt:38`
defaults `PLASIM_PLANET` to `p_earth`, and `exoplasim/scripts/build_model.py`
never sets it. A generic `p_exo.f90` sits beside it, unbuilt, and it is
materially better for this purpose: it drops `yplanet="Earth"` and `nplanet=3`,
zeroes the Earth Jan-1 `meananom0`, and adds `alv`, `als` and `tmelt` to
`planet_nl`. Verified from the run's own log rather than inferred: the diag
prints `YPLANET="Earth"`.

**II. A namelist "day" is a sidereal day, and this project feeds it Earth
days.** `day_24hr = 86400.0` (`plasimmod.f90:962`) is hardcoded, appears in no
namelist group anywhere in `plasim/src`, and cancels out of the
non-dimensionalisation: `dayseccheck` multiplies by it at `plasim.f90:1777-1780`
and `plasim.f90:1798-1810` divides by it again. The surviving convention is one
planetary rotation per namelist "day".

**III. Deleting Earth's input files does not yield a Vesper default.**
`configure()` correctly removes the shipped Earth boundary fields, and the diag
confirms codes 169, 210, 211, 709 and 903 all print `Init ... internally`. The
model does not then fall back to something neutral. It CONSTRUCTS an Earth
Arctic ice field, or reads an array nothing ever initialised.

**IV. Many of these constants have no namelist key at all.** A large share of
what follows sits in Fortran `parameter` statements or in module variables
absent from every namelist group. Those cannot be bracketed, swept, or declared
in `config/planet.yaml`, and moving one means a source edit and a rebuild of
every binary under rule 4.

---

# Live and load-bearing

## 1. Every cold start begins from a constructed Earth Arctic ice field and an ocean at its freezing point

`icemod.f90:119` leaves `xclsst` at its `-999` sentinel when code 169 is absent,
and `icemod.f90:230` then does `where (xclsst <= TFREEZE) xclicec = 1.0`, which
puts full ice cover on every ocean cell. `make_ice_thickness`
(`icemod.f90:1929-1982`) shapes it with `hmaxn=3.0`, `hmaxs=0.50`,
`cminn=0.1`, `cmins=0.25` and a fourteen-element monthly modulation table
credited in its own comment to "CCM3 report pp 127-129".

Separately, `oceanmod.f90:96` declares `real :: yclsst(NHOR,0:13)` with no
initialiser while its neighbours at lines 95, 97, 99 and 100 all have one.
`oceanmod.f90:285` does `yclsst = MAX(yclsst,TFREEZE)` and `:291-293` copies it
into the ocean temperature, so whatever was in memory is clipped up to the
freezing point.

Measured on 2026-08-24, first output bin of `run_2b20e3324bb0`:

    ts   ocean: min=214.832  max=271.900  mean=258.860
    sic  ocean: min=0.000    max=0.000    mean=0.000
    sit  ocean: NH mean=2.533   SH mean=0.139

The modelled ocean starts globally at `TFREEZE`. Every ocean cell carries ice
thickness in an 18:1 north-south ratio, which is Earth's Arctic-Antarctic
contrast imposed on a world with 32 degrees of obliquity and unrelated
geography, while ice COVER is identically zero: `xicecc` initialises to 0
(`icemod.f90:97`) and `mkicec` can only grow it. Thick ice carrying no albedo
is not a physical state, and neither is the field it replaced.

**What it costs.** The one long run escapes to 294.5 K with 0.07 per cent ice,
so the equilibrium survives. The exposure is `sweep_flux_earth`, which has not
been run: each arm will cold-start from a fully iced ocean, so the warm branch
is only ever approached from below and the hysteresis is never probed from the
ice-free side. An arm at 0.90 that fails to melt out would read as a snowball
bifurcation when it is an initial-condition artifact. Both fields are written to
the restart (`icemod.f90:995-997`, `oceanmod.f90:615`) and read back, so the
sentinel and the CCM3 thicknesses persist through every continuation. `yclsst`
is also emitted as ocean code 990, "clim. sst", which is meaningless.

## 2. The derived hyperdiffusion timescales are applied 25 per cent too weak

Mechanism II, in its most expensive instance. `check_consistency.py:728`
computes `base = math.pi*radius/(ntru*wind)/86400.0` and
`run_exoplasim.py:1031-1034` writes that straight into `TDISSD/Z/T/Q`. The
model reads it as sidereal days.

`dayseccheck` (`plasim.f90:1643-1656`) converts only when the values look like
days, `zmax < day_24hr/mtspd`; at `mtspd = 32` that threshold is 2700 and the
staged `TDISSZ = 1.1143` is far below it, so the conversion branch is taken.
Both routes to the error agree: reading day_24hr as cancelling gives a namelist
day of one rotation, and carrying it through explicitly gives a stored
coefficient 0.8 times the correct one.

The T21 rule gives tau_vorticity = pi * 7645464 / (21 * 5.94) = 192,552 s. The
model applies 240,678 s, exactly 108000/86400. All four timescales are 25 per
cent too long; the ratios between them survive, so the SHAPE of the damping is
right and only the level is wrong.

The confinement requirement still holds but with less margin than
`config/planet.yaml` states: solving `x = (dt/(2*kappa*tau))^(1/(gamma-1)) >=
0.6` needs gamma >= 14.78 at the nominal tau, which is the figure that block
quotes, and gamma >= 15.22 at the tau actually applied. `filter_power: 16`
clears both. So the integration is not in danger; what is wrong is that
`check_consistency.py` verifies a number the model does not use, and that
block's claim that the rule is "planet-agnostic by construction -- radius and
wind, no Earth days" is not true of what reaches the integration.

The same cancellation reaches `restim`, `tfrc`, `dampsp` and `taucool`. Of
those only `tfrc` is live; see finding 20.

## 3. The ocean's horizontal diffusion runs on Earth's radius

`oceanmod.f90:24`, `parameter(PLARAD=6.371E6)`, with the comment "Earth radius
(m)", used at `oceanmod.f90:1344` as `zfac=hdiffk(jlev)/plarad/plarad`.
`oceanmod` does `use resmod` rather than `use pumamod` (`oceanmod.f90:3`), so it
never sees the correct `plarad` that `planet_namelist` sets. It is a `parameter`,
so no namelist can move it.

Live and already exercised: `NHDIFF` and `HDIFFK` are written unconditionally by
`run_exoplasim.py:166-168`, and the CLIM-16 and CLIM-19 arms ran at 300, 1000
and 3000.

**What it costs.** (7.6452/6.371)^2 = 1.4396. Every arm applied 1.44 times the
tendency of its nominal diffusivity, so the bracket was really 432/1440/4318 and
the predicted rms heating table in
`exoplasim/notes/forcing-bundle-predictions.md` understates the model by the
same factor, 0.75/2.51/7.53 becoming 1.08/3.61/10.84.

This also resolves an open discrepancy. `predict_ocean_terms.py:99` builds its
operator with `a = radius_earth * EARTH_RADIUS_M`, which is Vesper's radius,
while claiming at line 58 to reproduce `hdiffo`'s discretisation, and its P2
eigenfunction check at line 80 validates the offline code against the analytic
Laplacian on that same radius -- so it structurally cannot catch this. The
prediction was +0.03 to +0.12 K and the measurement +0.13 to +0.19 K; dividing
the measurement by 1.44 lands it on the prediction. The 35 to 37 per cent
inter-realisation spread makes that a consistent explanation rather than a clean
reconciliation, but the factor is not in doubt.

## 4. A 5 m snow cap structurally neuters the glacier module that is switched on

`landmod.f90:55`, `dsmax = 5.00` m water equivalent, applied at
`landmod.f90:842-861`. The cap is applied to `dsnowz`, which is the same
accumulator `glaciermod.f90:355-358` builds ice orography from, so `glacieroro`
cannot exceed 5.0/0.85, about 5.9 m of ice. The module's stated purpose,
"capturing the effects of large ice sheets on atmospheric circulation"
(`glaciermod.f90:1-8`), is unreachable, and `where (dsnowz > 30.0) dglac = 1.0`
at `glaciermod.f90:211` is dead code.

5 m is PlaSim's Earth seasonal-snow ceiling, chosen for a model where perennial
ice arrives as a prescribed mask rather than being grown by this module.
`NGLACIER = 1` and `GLACELIM = 2.0` in the staged namelist, so the module is on.

**What it costs.** Everything above the cap is diagnosed as snowmelt
(`landmod.f90:860`) and pushed into `dwater`, so into soil water and runoff. Any
cell that would glaciate instead manufactures perpetual meltwater at up to the
full snowfall rate, inflating runoff and latent cooling and suppressing the
ice-albedo feedback exactly at the margin. This is the finding that most limits
what a run on this build can show, and it bears directly on
`notes/glacier-rough-pass.md`, whose areas are already an upper bound for other
reasons.

## 5. Cloud liquid water is distributed with an Earth scale height

`rainmod.f90:2005`, `zzh(:)=700.*ALOG(1.+dqvi(:))`, the CCM3 liquid-water scale
height, a length in metres fitted on Earth. It is used at `rainmod.f90:2008`
against mid-layer heights built at `:1990` and `:1995` as
`-dt*gascon/ga*ALOG(...)`, which correctly scale as 1/g. The heights shrink by
0.766 and the scale height does not. Both `700.` and the `0.00021` beside it are
bare literals with no namelist route.

Integrated over the T21/L10 sigma set at 25 kg/m2 precipitable water, the
Vesper-over-Earth liquid water per layer runs 7.2 at sigma 0.038, 4.1 at 0.12,
2.6 at 0.21, 1.9 at 0.32, 1.5 at 0.44, 1.2 at 0.57, 1.0 at 0.70, 0.89 at 0.82,
0.81 at 0.92 and 0.78 at 0.98.

**What it costs.** Column-integrated liquid water barely moves, 5 to 7 per cent
low, but mid and upper cloud carries 1.5 to 4 times the calibrated water at the
same cloud fraction while the lowest three layers carry 11 to 22 per cent less.
That is the wrong direction on both terms that matter: thicker high cloud raises
the longwave trap, thinner low cloud lowers the shortwave reflection. `dql` is
not diagnostic; it sets shortwave cloud optical depth (`radmod.f90:2465`) and
longwave cloud emissivity (`radmod.f90:3251`, `:3270`). The consistent
correction is to scale the 700 by g_earth/g, to about 536, since liquid water
should track vapour and vapour's geometric scale height is 0.766 here.

## 6. Land longwave surface emissivity is hardcoded to exactly 1.0

`radmod.f90:3237`, `zeps(:) = dls(:) + 0.98*(1.-dls(:))`. The modelled ocean
gets 0.98 and the modelled land gets a perfect blackbody. No namelist key
exists.

This is a generic Earth-GCM simplification and it costs more here than on Earth.
`config/planet.yaml` names evaporite and playa clastics as barren classes
covering about a quarter of this land, and salt crust and quartz-rich clastics
run 0.90 to 0.95 broadband in the 8 to 14 micron window.

**What it costs.** At eps = 0.92 and a modelled surface at 300 K the surface
emits 37 W/m2 less and reflects about 26 W/m2 of downwelling that the code
discards outright, since `zeps = (1-zeps)*dftd(:,NLEP)` is identically zero over
land. Net about 11 W/m2 over that terrain, which scaled by a quarter of the land
fraction is **0.4 to 1.1 W/m2 planetary mean**, bracketed on the emissivity
range. That is larger than the adiabatic term `energy_fixer` was adopted to
compensate.

## 7. The snow-over-forest albedo is one spectrum-blind factor applied to both bands

`landmod.f90:229-236`:

    albsmaxf1 = 0.5*albsmax1 ; albsmaxf2 = 0.5*albsmax2
    albsminf1 = 0.75*albsmaxf1 ; albsminf2 = 0.75*albsmaxf2

used at `landmod.f90:544-561` to set `dsalb(1,:)` and `dsalb(2,:)`, which are
what the radiation reads at `NSTARTEMP=1` (`radmod.f90:2951-2952`).

The endmembers `albsmax1/2` are derived from the k25v spectrum and are correct.
The canopy masking on top of them is not: 0.5 and 0.75 are Earth-Sun broadband
ratios, from upstream's `albsmaxf=0.4` over `albsmax=0.8`, applied identically
either side of 0.75 micron. This is the class
`inherited-earth-constants.md` finding 1 named for codes 175 and 176, appearing
here in the SNOW albedo, where nothing has corrected it.

Measured on 2026-08-24 from the staged `orogen_T21_surf_0212.sra`: 2048 cells,
min 0.0000, max 0.4988, nonzero on 997 cells with a mean of 0.3120 over those.
The field is live, not the uniform-land default.

**What it costs.** With `albsmax1 = 0.96819`, `albsmax2 = 0.58529` and band
weights 0.38237/0.61763, anchoring on band 1 with this project's own vegetation
band pair implies a canopy cover of 0.542, which gives a band-2 forested
endmember of 0.390 against the code's 0.293, so 0.097 too dark at full canopy.
Anchor on band 2 instead and band 1 comes out about 0.16 too bright; one factor
cannot serve both. At the land-mean forest fraction that is about 0.019 too dark
broadband and 0.030 at the maximum, one-signed warm, and the band carrying the
larger half of the error is the one this star puts 0.618 of its surface
shortwave into.

## 8. The critical relative humidity is tuned to an Earth-sized grid box

`rainmod.f90:92`, `rcrit(:)=MAX(0.85,MAX(sigma(:),1.-sigma(:)))`, the threshold
for all non-convective cloud at `rainmod.f90:1971` and `:1977`. It is a subgrid
humidity-variance parameter, so it is a function of grid box AREA, and 0.85 is a
T21-on-Earth tuning. This planet's radius is 1.20 times Earth's, so a T21 box is
1.44 times the area and carries more subgrid variance, which means `rcrit`
should be LOWER.

`RCRIT`, `RCRITMOD` and `RCRITSLOPE` are namelist keys (`rainmod.f90:79-81`) and
this project sets none of them.

**What it costs.** The model produces systematically too little stratiform cloud
at a given grid-mean relative humidity, on the single largest lever this scheme
has over the planet's albedo.

## 9. Precipitation re-evaporation is calibrated to an Earth column and an Earth timestep

`rainmod.f90:31`, `gamma = 0.01`, used at `rainmod.f90:2238`, `:2249`, `:2260`
and `:2271`. `nevapprec = 1` by default (`rainmod.f90:20`), so it is on.

Two defects in one constant. The physical content is the fraction of a layer's
sub-saturation that falling precipitation fills in transit, which is
proportional to residence time dz/v_t. Here dz goes as 1/g and is 0.766 times
Earth's, and terminal velocity goes as sqrt(g) at the same air density and is
1.143 times, so residence time is **0.67 times Earth's** and the Earth tuning
over-evaporates by about 1.5 times in every layer. Separately, the evaporated
AMOUNT per step carries no `deltsec2`, so the RATE is inversely proportional to
the timestep.

**What it costs.** The first half lands on P minus E over land, which is the
numerator of the carve criterion. The second means this term's strength varies
by a factor of two across `resolution_timestep_minutes`, from 45 down to 22.5,
with no physical content -- so any comparison across the resolution ladder in
loop D is comparing two different re-evaporation strengths.

## 10. Salinity reaches the model through four compiled constants, not one

`config/planet.yaml` states that salinity "reaches the model through one number,
the freezing point of sea water". Four Fortran `parameter`s carry it and no
namelist key exists for three of them.

| constant | site | what it is | across S = 20 to 40 |
| --- | --- | --- | --- |
| `TFREEZE` | `icemod.f90:54` | freezing point | handled; -0.04 to -0.36 K |
| `CRHOS = 1030.` | `icemod.f90:21`, `oceanmod.f90:20` | density "at S=34.7" | 1016 to 1032 kg/m3 |
| `CPS = 4180.` | `oceanmod.f90:22` | labelled sea water; this is FRESH water's value | 4070 to 3950 |
| `CLFI = 3.28E5` | `icemod.f90:30` | fusion, depressed by Earth brine content | toward 3.337e5 as it freshens |

**What it costs, in two separable ways.** The mixed-layer heat capacity
`CRHOS*CPS` used at `oceanmod.f90:401` and `:816` is 4.305e6 J/m3/K against sea
water's real 4.10e6 at the declared salinity, so it is 5.1 per cent too large
and the DECLARED `mixed_layer_depth_m: 50.0` is thermally 52.5 m of sea water.
`scripts/error_budget.py` carries that depth as a structural item stated in
metres, and the number it carries is 5 per cent off from what it means.

Second, the snow-ice flooding threshold at `icemod.f90:1349-1352` is
`1.E3*xsnow > (CRHOS-CRHOI)*xiced`, a DIFFERENCE of 1030 - 920 = 110, so a 1 per
cent density error is a 10 per cent threshold error. Across the declared bracket
the true margin runs 96 to 112, so -13 to +2 per cent about the compiled value: a
fresher ocean floods and converts snow to ice substantially earlier, trading
insulating snow for conducting ice and accelerating basal growth. That is a
sea-ice feedback CLIM-17's bracket is meant to span and currently cannot.

## 11. The ozone profile is placed at a fixed geometric height

`radmod.f90:81-82`, `bo3 = 20000.` and `co3 = 5000.`, both in metres, applied at
`radmod.f90:1987-1996`. `mko3` builds its height coordinate hypsometrically with
the model's own `ga`, correctly, and then places the profile against a fixed
geometric 20 km. This world's scale height is 0.766 of Earth's, so 20 km is 3.6
scale heights here against 2.7 on Earth.

The modelled photochemical ozone maximum is set by pressure-like conditions --
ultraviolet optical depth and three-body recombination density -- not by
geometric altitude, so holding the height fixed is the wrong invariant. Holding
pressure fixed means `bo3` about 15300 and `co3` about 3830.

The keys exist and are unreachable as the tree stands:
`vendor/exoplasim/exoplasim/__init__.py:2927-2935` accepts `ozone` as a dict
with `height` and `spread` and writes `BO3`/`CO3`, while
`run_exoplasim.py:1894` passes `ozone=bool(atmosphere["ozone"])`, forcing the
boolean path. Both files need a line.

**What it costs.** Reproducing the L10 sigma levels and the `mko3` recursion
with a 230 to 288 K profile, the fraction of the modelled ozone column deposited
in the top model layer goes from 0.616 at Earth gravity to **0.794** here, and
layer 2's share halves from 0.202 to 0.105. This redistributes shortwave heating
rather than changing the column, so the budget effect is small and consistent
with `exoplasim/notes/ozone.md` calling it second-order for the surface; the
top-layer heating rate rises about 29 per cent, and that is the layer setting
the modelled tropopause temperature and static stability. Distinct from the
`a1o3`/`aco3`/`toffo3` shape issue that note already records: this one is about
gravity.

## 12. The sea-ice albedo ramp keeps an Earth broadband slope under star-derived endpoints

`seamod.f90:156-158` and `:270-272`:

    dsalb(b,:) = doceanalb(b)*(1-dicec) + dicec*AMIN1(dicealbmx(b), dicealbmn(b) + 0.025*(273.-dts))

`dicealbmn` and `dicealbmx` are star-weighted per band by `radmod.f90:725-756`
and independently verified by `analysis/ice_albedo.py`. That work stopped at the
endpoints. The slope 0.025 per K is Earth's, tuned to Earth's broadband span
`albice - 0.5 = 0.25`, which is exactly 10 K, and neither it nor the 273.0
anchor is a namelist key. `landmod.f90:436-441` does the same job correctly with
a FRACTIONAL ramp, `(albgmax-albgmin)*(dts-263.16)/(tmelt-263.16)`, which is the
structure this should have.

**What it costs.** Band 1 spans 0.6203 to 0.8702, a span of 0.250, and saturates
10.0 K below freezing; band 2 spans 0.3395 to 0.4729, a span of 0.133, and
saturates 5.3 K below. The two bands therefore reach maximum albedo 4.7 K apart,
a gap set entirely by how far the re-weighting stretched each endpoint and by
nothing about the modelled ice. At 268 K band 2 reads 0.4645 against a
consistently scaled 0.4062, so +0.058 in band 2 and about +0.036 broadband. Near
zero in the global mean, because most modelled ice is colder than either
saturation point, and several W/m2 locally in the 263 to 273 K band where the
ice-albedo feedback lives. One tell that 273.0 is a magic number rather than a
physical anchor: it is neither `TMELT = 273.16` nor `TFREEZE`, so at the melting
point the formula returns 0.004 BELOW its own declared minimum.

## 13. The free-convection ocean transfer coefficient embeds Earth's gravity

`fluxmod.f90:254-255`, the Miller, Beljaars and Palmer (1992) unstable-over-ocean
branch:

    zrifh(jhor)=(1.+(0.0016*zdth**(1./3.)/SQRT(zabsu2(jhor))/zkblnz2)**1.25)**0.8

In the free-convection limit this collapses to a transfer velocity of
`0.0016 * dtheta**(1/3)`, the classical four-thirds flux law, in which the
0.0016 IS the `(g/theta)^(1/3)` group evaluated at 9.81. This world needs
`0.0016*(12.81/9.81)^(1/3) = 0.00175`, **9.3 per cent higher**. Unlike the
Charnock case there is no compensating length scale, because the four-thirds law
carries no boundary-layer depth. No namelist key exists.

**What it costs.** `dtransh` feeds both `mkshfl` (`fluxmod.f90:534`) and
`mkevap` (`:664`), so this is modelled ocean sensible heat and evaporation
together. At the lowest level height of about 111 m and a sea roughness of order
3e-5, the neutral coefficient is about 7e-4, and at 2 m/s with a 1 K air-sea
difference the free-convection term already nearly doubles the transfer
coefficient; the sensitivity there is about 0.54 and rises to 1.0 in the fully
calm limit. So a 5 to 9 per cent low bias in latent and sensible heat over the
calm tropics.

Structural and worth recording beside it: there is no free-convection branch
over land at all (`fluxmod.f90:257`), so calm daytime fluxes over the playa and
evaporite quarter get no gustiness enhancement.

## 14. One Earth-average soil thermal pair, and the conductivity cannot be set

`landmod.f90:73`, `soildiff = 1.8` W/m/K, and `landmod.f90:76`,
`soilcap = 2.4E6` J/m3/K, applied uniformly at `landmod.f90:697-698`. These are
moist mineral soil, Earth's global average, on a project that already derives
albedo and roughness per cell from lithology.

`SOILCAP` is a namelist scalar (`landmod.f90:172`), reachable through `cpsoil`
and not passed. **`soildiff` is not in the namelist at all**, and neither are
`snowdiff`, `sicediff`, `snowcap`, `sicecap` or `rhosnow`. Neither can be a
field; both are scalars.

**What it costs.** Thermal inertia I = sqrt(k*rho*c) gives 2078 J/m2/K/s^0.5 for
the default against roughly 625 for a dry playa or salt crust, so about a quarter
of the modelled land carries **3.3 times too much thermal inertia** and its
diurnal and seasonal surface temperature swing is damped by roughly that factor.
That feeds evaporation, which feeds P minus E, and it feeds the dust emission
threshold.

---

# Live, bounded or diagnostic

## 15. The archived surface albedo is the Earth-Sun broadband value, not the one the radiation used

`landmod.f90:28-33` declares non-spectral `albsmin=0.4`, `albsmax=0.8`,
`albsminf=0.3`, `albsmaxf=0.4`, `albgmin=0.6`, `albgmax=0.8`. `landini`
re-derives only the banded copies from the spectrum (`landmod.f90:225-236`); the
six scalars are never touched. `landstep` computes `dalb(jhor)` from them
(`:544-557`, `:601-603`), and since `radstep` runs before `surfstep` in the
timestep (`plasim.f90:4048` then `:4082`), radmod's correct overwrite of `dalb`
at `radmod.f90:2943` is discarded by landmod every step. The landmod value is
what `outmod.f90:2641` accumulates and `:695` writes as code 175, which is in
both `REGULAR_CODES` and `SNAPSHOT_CODES`.

Radiation itself reads `dsalb` and is unaffected, so this is archival only. The
stored field reads +0.068 over cold snow, +0.057 over snowy forest and +0.061 to
+0.068 over cold and melting glacier. Any future surface energy closure or map
that reads 175 will be wrong by that much, one-signed, only over the cryosphere.
Code 175 also collides with `dalbcl1` on input (`surfmod.f90:252`), so the same
code means two different things in and out.

## 16. Mean sea-level pressure is reduced with Earth's standard-atmosphere lapse rate

`vendor/exoplasim/exoplasim/pyburn.py:217`, `RLAPSE = 0.0065`, used at `:2001`
and `:2004` and mirrored at `:2841` and `:2844`, with the ECMWF cold-surface
correction `tstar[tstar<255.0] = 0.5*(255+tstar[tstar<255.0])` at `:2003`
beside it. Everything else in that block, `gascon` and `gravity`, is planet-aware
and passed down from `configure()`; `RLAPSE` is the sole exception in the whole
postprocessor, and `pyburn.postprocess` has no lapse-rate argument.

This world's dry adiabat is 12.75 K/km, so the same fraction of it that 6.5 K/km
is of Earth's gives about 8.5 K/km. The reduction exponent `alpha = gascon*RLAPSE/gravity`
goes from a consistent 0.1904 to the coded 0.1456, which is +5.3 hPa at 2 km of
terrain and +10.3 hPa at 4 km, one-signed and orography-shaped.

Code 151 is in `REGULAR_CODES` and `SNAPSHOT_CODES`, so `psl` is written to
every climatology, and nothing in the project reads it. A trap rather than a
current error: the field is written, carries a CF standard name, and looks
authoritative. Same class as `inherited-earth-constants.md` finding 3, which
found 6.5 K/km in three project scripts and not in the postprocessor.

## 17. Earth's lapse rate and a 288 K surface set every cold start's initial profile

`setzt` (`plasim.f90:2202-2259`) builds the restoration and initial temperature
profile from `tgr = 288.0` (`plasimmod.f90:314`), `ALR = 0.0065`
(`p_earth.f90:42`) and `dtrop = 12000.0` (`plasimmod.f90:312`). All three are
namelist-reachable, through `plasim_nl` for the first and third and `planet_nl`
for `ALR`, and this project sets none of them. `setzt` is reached through
`initfd`, which `plasim.f90:487-491` gates on `nrestart == 0`, so this is the
cold-start path only.

`restim = 0` so `damp = 0` and `sr` is an initial state rather than a
restoration; `dtep = dtns = 0` so the latitudinal structure is null. This
world's dry adiabat is 12.76 K/km against Earth's 9.76 and its scale height is
0.77 times, so a 6.5 K/km profile with a 12 km tropopause is a poorer first
guess here than it is on Earth. It washes out, so this costs spin-up rather than
an answer -- but `config/planet.yaml` records that the resolution-ladder
stability boundary was measured entirely on cold starts, and that boundary is
the ceiling for exactly this transient.

Second and smaller use: `plasim.f90:569` does
`zdlnp = zmeanoro/gascon/tgr; psurf = exp(log(psurf)-zdlnp)`, so the declared 1
bar is reduced for mean orography using Earth's 288 K. At a mean orographic
geopotential of order 2500 m2/s2 that is about 200 Pa, 0.2 per cent of surface
pressure, but it means the realised global-mean surface pressure depends on an
undeclared Earth number.

## 18. The model's calendar is a 360-day Earth year that `calini` never updates

`calmod.f90:19-22` declares `n_days_per_month=30`, `n_days_per_year=360`,
`m_days_per_year=360`, `m_days_per_month=30`. `calini` (`calmod.f90:44-99`)
copies `n_days_per_month`, `n_days_per_year`, `n_start_step`, `ntspd`,
`solar_day` and `mpstep` -- and NOT `m_days_per_year` or `m_days_per_month`,
which `plasim.f90:1506-1507` has just set to 183 and 15.

So `mtspd` in calmod is `146*40/360`, truncating to 16, and the calendar year is
`360*16 = 5760` steps against the orbit's `n_steps_per_year = 5850`. The
calendar drifts 1.54 per cent per orbit against the season and is a full year
out after about 65 orbits. `tcalday` becomes 43200 s, so a "calendar day" is
half a 24-hour day and `step2cal30` reports hours 0 to 11. `cal2step` and
`step2cal30` are not inverses: `plasim.f90:1585` encodes with `mtspd = 32`, 146
and 15, and `step2cal30` decodes with 16, 360 and 30, so year 1 month 1 day 1
encodes to step 4672 and decodes as `23-Oct-0001`. `outmod.f90:119` and `:147`
then pack `ihead(3)` from that while `ihead(8)` reports pumamod's
`m_days_per_year = 183`, so one eight-word header carries two different days per
year.

**Bounded today.** `pyburn.py:492` and `:539` take the netCDF time axis from
`header[6]`, the raw step counter, and radiation takes its orbital phase from
`n_steps_per_year` independently (`radmod.f90:1848`, `:1979`), so nothing this
project reads goes through the broken calendar. It becomes live the moment
anyone enables climatological ozone (`radmod.f90:1999`), t-nudging or flux
correction (`miscmod.f90:264`, `:321`), the land surface annual cycle
(`landmod.f90:1490`, `:1519`), or prescribed ice and SST (`icemod.f90:1883-1919`,
`oceanmod.f90:745`, `:765`) -- all of which interpolate through `step2cal30`.

Related trap in the same machinery: `n_days_per_year == 365` silently switches
the whole model to Earth's Gregorian calendar, with `mondays`, the 400/100/4 leap
rule and `Jan..Dec` names (`calmod.f90:107`, `:138`, `:277`, `:387`, `:463`).
`N_DAYS_PER_YEAR` is derived from the flux and the rotation period at
`run_exoplasim.py:238`, so nothing structurally prevents a future world landing
on 365, and there is no guard and no message.

## 19. The glacier persistence test is one run segment, which here is half an Earth year

`glaciermod.f90:42` initialises `persistflag(NHOR) = .TRUE.`; `glacierstep`
(`:588`) only ever clears it; `glacierstop` (`:614`) converts whatever survives
to `dglac = 1.0`. **`persistflag` is not written to the restart** --
`glacierstop` writes only `groundsg`, `dglacsg`, `doro` and `dglac` at
`:622-628`.

The module header says the criterion is "some snow cover continuously for an
entire year" (`glaciermod.f90:10-11`). What it is, is continuously for one
invocation of the model. This project runs one orbit per invocation, so the
persistence evidence a modelled glacier needs is half what the module was
written to require, and it resets to a clean slate at every segment boundary.
The threshold `GLACELIM` is reachable and is set to 2.0, which is ExoPlaSim's
Earth snowline default restated; the DURATION is not reachable at all, because
it is the segment length.

## 20. An undeclared Rayleigh sponge of 20 and 100 days on the top two layers

`plasim.f90:1437-1439`, run before the namelist read:

    tfrc(1) = 20.0*day_24hr*frcmod
    tfrc(2) = 100.0*day_24hr*frcmod

`TFRC` and `FRCMOD` are in `plasim_nl` and this project sets neither. Note the
asymmetry: the NLEV=20 equivalent (`plasim.f90:1459-1469`) requires
`nrdrag = 1` to switch on, while the NLEV=10 case, which is the configuration
this project runs, is unconditional. By mechanism II these are 20 and 100
sidereal days, so 25 and 125 Earth days.

An inherited Earth-tuned top-of-model sponge, live in every run, on the layers
at sigma 0.025 and 0.094 -- 2.5 and 9.4 kPa on a 1 bar surface -- appearing in
none of this project's declared parameters, unlike the filter, the
hyperdiffusion and the timestep, which are all argued in `config/planet.yaml`.
Whether 20 days suits a 30-hour rotator at 1.31 g is the kind of question this
project's conventions say is answered rather than inherited.

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

`fluxmod.f90:26`, `vdiff_lamm = 160.`, the asymptotic Blackadar mixing length in
metres, used at `:783` and `:805-806` with `ztscal = 250.` hardcoded at `:762`:

    zzlev = -gascon*ztscal*ALOG(sigmah(jlev))/ga
    zmixm = vdiff_lamm*vonkarman*zzlev/(vdiff_lamm+vonkarman*zzlev)

`zzlev` scales as 1/g and `vdiff_lamm` does not, so the surface-layer to
free-atmosphere crossover sits at a fixed geometric 400 m rather than at a fixed
fraction of a boundary layer that is 0.77 times as deep here. Diffusivity goes
as the square of the length: at the lowest interface, sigmah about 0.966, the
length ratio is 0.845 so K is 71 per cent of the Earth-equivalent, and by
sigmah 0.7 the deficit has fallen to 8 per cent. Under-mixing concentrated
exactly where surface flux is redistributed. `VDIFF_LAMM` is a namelist key
nothing sets.

## 22. Snow masks the surface at 0.01 m water equivalent regardless of what it buries

`landmod.f90:406-410` and `:556-561`,
`dalb = dalbclim + (zalbsnow - dalbclim)*dsnow/(dsnow+0.01)`. Half the snow
albedo effect at 0.01 m water equivalent, about 3 cm of modelled snow at
`rhosnow = 330`. This is a snow-cover-fraction parameterisation whose scale
should track the roughness of what is being buried, and this project's own
staged `orogen_T21_surf_0173.sra` gives land roughness from 0.0038 m to 7.57 m
with a land mean of 2.0 m.

A single 3 cm masking depth turns a 7.6 m rough cell half-white on the first
snowfall. One-signed, too reflective too early, and it applies to the spectral
`dsalb` as well as to `dalb`. The constant is inline in two places with no
namelist route.

## 23. Convective cloud cover is an Earth mm/day regression evaluated in 30-hour days

`rainmod.f90:1926-1928` with `zcca=0.245`, `zccb=0.125` at `:1905`:

    zrfac = solar_day * 1000.0 ! convert m/s into mm/day
    zcctot(:)=zcca+zccb*log(dprc(:)*zrfac)

Credit where due: the fork replaced a hardcoded 86400 with `solar_day` here, and
the three uncompiled convection variants still carry `zrfac=86400.*1000.`
(`rainmod_bm.f90:1246`, `rainmod_mca.f90:837`, `rainmod_kuo_old.f90:1146`). But
the regression was fitted on Earth precipitation rates in Earth days, so neither
day length is defensible on the fit's own terms.

Using this world's solar day inflates the log argument by 1.26 and adds a
uniform **+0.029 to convective cloud fraction** wherever the result is not
clamped to [0.05, 0.8]. The point is not which choice is right but that a free
choice here is worth 3 per cent absolute cloud cover.

## 24. There is no saturation-over-ice branch, and the remedy is blocked by mechanism I

`p_earth.f90:44-46` sets `ra1=610.78`, `ra2=17.2693882`, `ra4=35.86`, the
Magnus-Teten coefficients over LIQUID water, and that single formula is used at
every `zqsat` site in `rainmod.f90` and at `seamod.f90:254`. Meanwhile the
latent heat correctly switches to `ALS` below `TMELT` at `rainmod.f90:767`,
`:886`, `:980`, `:1083` and at `fluxmod.f90:693-700`, and the Clausius-Clapeyron
derivative is the liquid one multiplied by L_s/cp.

Saturation over ice is about 25 per cent below saturation over liquid at 250 K,
so cold-cloud condensation is systematically over-produced and the
thermodynamics is internally inconsistent with its own latent heat.

The ice coefficient set that the Magnus formula needs over ice is
610.66/21.875/7.65, measured on Earth's water. `alv`, `als` and `tmelt` are now
in `planet_nl` on the module that compiles, so the latent heats and the
switching temperature are settable, but `ra1`, `ra2` and `ra4` are one set and
the fix needs a second BRANCH rather than a different single set. It is a source
change either way. world-ako.

## 25. Snow and glacier densities are Earth compaction values

`landmod.f90:72`, `rhosnow = 330` kg/m3, and `glaciermod.f90:47`,
`rhoglac = 850` kg/m3. `rhosnow` converts water equivalent to physical snow
thickness for the top-layer heat capacity blend (`landmod.f90:734`, `:741`,
`:784`, `:791`), so it decides how much the modelled snowpack insulates the
soil; `rhoglac` converts water equivalent to ice thickness for the orography
(`glaciermod.f90:355`). Both are settled densities set by overburden compaction,
which scales with gravity, and this world's is 1.306 times Earth's, so settled
snow would be denser and firn would densify faster.

The fork got the GRAVITY right at `glaciermod.f90:344-351`, with a comment
explaining why, and then multiplies it by an Earth-compacted density. Worth
roughly 10 to 20 per cent on snow insulation thickness and on glacier relief.
Neither is in any namelist.

## 26. The lead-closing scale and the thickness clamp are Earth Arctic tunings

`icemod.f90:1152`, `zh0 = 0.5` in `mkicec`, credited to Hippler 1979 at
`:1163`, with the same 0.5 m in `mkicecf` at `:1205`. It is the whole of the
model's lead parameterisation: a cell's compactness closes with an e-folding of
0.5 m of new ice growth, and `icemod.f90:746-750` then thresholds it at
`thicec = 0.5` into a hard mask. So `zh0` alone sets how much modelled ice must
grow before a cell flips to iced for albedo and roughness, roughly 0.35 m. On a
world whose year is half Earth's the ice grown per season differs, and this
Arctic-calibrated growth scale is what converts it into cover. `thicec` is a
namelist key and unset; `zh0` is not a key at all.

Beside it, `xmaxd = 9.0` m (`icemod.f90:58`) is Earth's Arctic multi-year
maximum, and its enforcement is worse than a clamp: `icemod.f90:705-718` melts
the excess and `getiflx` (`:1988-2020`) takes the GLOBAL area-weighted sum of
that heat and redistributes it onto every other cell with ice below 9 m, which
is a non-local heat transport with no physical carrier. Dormant at the current
climate -- thickness peaked at 4.46 m in spin-up and 0.78 m at equilibrium, so
the branch never fires -- and live in any cold sweep arm.

## 27. Runoff velocity constants are Earth-fitted, with gravity entering through the geopotential

`landmod.f90:1185-1186` and the identical pair at `glaciermod.f90:303-304`,
`zcvel = 4.2` and `zcexp = 0.18`, used at `landmod.f90:1328-1347` and
`glaciermod.f90:541-561` as `u = zcvel/zdx * |dh/dx|^0.18`. `zoro` there is
GEOPOTENTIAL, not elevation, so the slope term carries `ga`: the same
topographic slope gives 1.306 times the value here, worth 1.306^0.18 = 1.049 on
the velocity. Separately, open-channel velocity scales as sqrt(g), so a
physically consistent `zcvel` would be about 1.14 times larger. The pit-filling
increment `zoron = 1. + MIN(...)` is likewise 1 m2/s2, so 0.078 m of elevation
here against 0.102 m on Earth.

`oroini` runs after `roffini` in `surfini` (`surfmod.f90:479-480`), so the
glaciermod copy is the one that survives. A few per cent on river transit,
recorded because it is the only routing timescale in the model and it is
unreachable.

## 28. Every run's log identifies the planet as Earth

`p_earth.f90:88-105`. `print_planet` has no namelist and nothing parses it, so
this is diagnostic only -- but it goes into `plasim_diag` beside numbers that are
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
temperature, perihelion and aphelion are printed directly above this world's
correct mean radius, gravity and irradiance. Only `plarad`, `ga`, `gsol0` and
`sidereal_day` in that table are the real values. This world is 1.881 Earth
masses at 1.20 Earth radii, so a density near 7600 kg/m3 against the 5520
printed, and an equilibrium temperature near 245 K at A = 0.3 against 247.3 that
agrees by coincidence. Not merely stale: a different planet, in the artifact a
reader consults to confirm what ran.

---

# Dormant, and what each waits on

None of these affects a run today. Each is recorded because it fires on a change
this project already intends to make, and because a dormant Earth constant reads
exactly like a live one to whoever throws the switch.

**Dust would settle about 1790 times too slowly.** `aeromod.f90:43-44` declares
`apart = 50e-9` m and `rhop = 1000` kg/m3, a photochemical haze grain at water
density, and both are in `aero_nl`. `run_exoplasim.py:1154-1252` writes every
key in `prov["namelist_values"]`, and
`aeolian/scripts/build_dust_source_fields.py:98-129` builds that dict with the
fourteen `DUST*` keys and NEITHER `APART` NOR `RHOP`. The correct values exist
unread in `exoplasim/data/dust/vesper_dust_aerosol.provenance.json` as
`apart_m = 2.2068e-06` and `rhop_kg_m3 = 2600.0`, derived for this planet's
gravity by `dust_aerofile.py:158`. Through `vels` at `aerocore.f90:1185` the
radius contributes 1948, the density 2.602 and the Cunningham factor 0.354, so
1793 net; sedimentation is the only removal term active by default, so the
burden falls back on the timestep-dependent 99-per-cent-per-step scrub at
`aerocore.f90:944`. Fires on `model.dust_emission`. Note that
`aeolian/notes/in-model-dust.md:528` asserts the sidecar carries `APART` and
`RHOP` in `namelist_values`, and the code does not do this. A new instance of
the class in `aerosol-particle-radius.md`, not the audited one: that audit is
about radmod's copy diverging from aeromod's, and this is about neither copy
ever being set.

**SIMBA is Earth-fitted throughout, with two latent traps.** Unreachable at
`NVEG = 0` and one namelist key away. `rlue = 3.4E-10` kg C/J
(`landmod.f90:63`) is applied to broadband downward shortwave (`simba.f90:387`)
rather than to the photosynthetically active fraction, which under this star is
roughly 0.30 to 0.33 against about 0.45 for the Sun, so light-limited modelled
GPP would be overstated by about 1.4 times; the code says so itself at
`simba.f90:456-459`. `ct_crit = 5.0` gives a linear ramp from 273.16 K and then
flat forever, with no temperature optimum and no high-temperature cutoff;
`q10 = 2.0` is referenced to 283.16 K; `co2_ref = 360` ppm and `co2_sens = 0.3`
are Harvey (1989)'s Earth beta factor; `zlaimax`, `cveg_k`, `cveg_l`, `cveg_a`
and `cveg_f` are Earth biome fits; `valb_min`, `valb_max` and `vsalb_min` are
Earth-Sun broadband albedos, and `vsalb_min` is a THIRD independent forest-snow
albedo constant inconsistent with both finding 7 and anything derived from this
star. `tau_veg = 10` and `tau_soil = 42` years are scaled correctly and
planet-aware at `landmod.f90:296-297`, which is the trap in the other direction:
it converts Earth-calibrated turnover times into VESPER years, so 10 and 42
become 5 and 21 Earth years of physical turnover. The two latent traps:
`simba.f90:484` overwrites `dwmax` every timestep from SIMBA's own Earth range,
silently discarding the pedology-derived code-229 field the whole
bootstrap-to-baseline sequence exists to produce; and `simba.f90:486` sets
`dalb` but not `dsalb`, so under `two_band_albedo: true` the radiation reads
`dsalb` only and SIMBA's vegetation albedo would have no radiative effect at all
while its roughness and bucket would. `nveg=2` and `twobandalbedo=True` are
incompatible and nothing says so.

**Flux correction relaxes toward the constructed Earth ice field.** `nfluko` is
the model's actual q-flux mechanism. `nfluko=1` reads codes 709 and 903, which
do not exist here so they stay zero. `nfluko=2` (`mkflukoi`, `icemod.f90:1582`)
relaxes ice toward `xcliced2`, which is finding 1's CCM3 field, and `addfci`
(`:1653-1654`) compares SST against the interpolation of the `-999` `xclsst`.
`config/planet.yaml` explicitly wants a q-flux bracket and
`analysis/error_budget.json` carries no-q-flux as a structural term, so this is
directly on the path of work already planned.

**The resolution-tuned radiation constants are pinned at T21's row by two
independent mechanisms.** `radmod.f90:945-1005` sets `tswr1`, `tswr2`, `tswr3`
and `th2oc` per (NTRU, NLEV), and every branch is guarded by
`if(NDCYCLE==1) jtune=0`. `ndcycle` defaults to 1 (`radmod.f90:196`) and nothing
sets it, so every branch of that table is dead and the module defaults stand at
every rung -- and those defaults ARE the T21/L10 row. Second mechanism:
`run_exoplasim.py:171-175` writes `TSWR3 = 0.0055 * cloud_absorption_scale`,
hardcoding T21's base. Correct today at T21. At T42 the model would use `tswr1`
0.077 against the table's 0.089, `tswr3` 0.0055 against 0.0048, and `th2oc`
0.024 against 0.0285, so 16 per cent on the longwave water-vapour continuum.
`th2oc` is a new member of the cloud-tuning class of
`inherited-earth-constants.md` finding 2 and is not on the corrected list.

**Berger's Milankovitch series and an Earth day-80.5 equinox.**
`radmod.f90:3773-4177` computes Earth's orbital elements from a year AD, reached
from `radini` when `nfixorb == 0`; `run_exoplasim.py:1885` sets
`fixedorbit=True`, so it is unreachable. `radmod.f90:4295`,
`lambm = lambm0 + (calday - ve/365.)*2.*pie` with `ve = 80.5`, is Earth's
Jan-1-to-equinox phase and does not scale with the model calendar; dormant
because `keplerian=True` routes `solang` to `gen_orb_decl`, which takes its
phase from `mvelpp` and `meananom0r` and is correct and config-reachable. Both
become live on one keyword.

**Betts-Miller's Earth adjustment timescales.** `rainmod_bm.f90:414-415`,
`ztaud=7200.` and `ztaus=14400.`, the canonical two-hour deep and four-hour
shallow Earth relaxation. Fires on swapping `RAINMOD` in the build.

**Hurricane diagnostics are Earth-empirical end to end.** Dormant at
`nstormdiag = 0`. `hurricanemod.f90:62-74` is a wall of Earth tropical-cyclone
thresholds. Two are specifically broken off Earth rather than merely unfitted:
`LAVTHRESH=1.2e-5` per s absolute vorticity, against a planetary vorticity 0.8
times Earth's at 30 hours, so 25 per cent too strict; and `SIZETHRESH = 30`
cells against a T21 cell 1.44 times Earth's area, so a 44 per cent larger storm.
The GPI formula at `:1355-1358` carries its Earth normalisations inline. `RD` is
correctly taken from `gascon` at `:181`, while `CPD=1005.7` at `:52` is
hardcoded and not tied to `acpd`. The positive form: this diagnostic is not
meaningful on this world and should stay off; if it is ever wanted, the
thresholds are a re-derivation and not a namelist tweak.

**LSG and the coupler are Earth-hardwired and not compiled.**
`plasim/CMakeLists.txt:182-183` links `src/lsgmod.f90`, a four-line stub, and
`cpl_stub.f90`. `lsg/src/lsgmod.f90:3176-3199` and `:3263-3265` carry
`g=9.80665`, `erdrad=6371000.00`, `erdrot=four*pi/86164.`, `rhonul=1030.`,
`tfreez=-1.91`, `cp=4180.` and `entmel=80.*cp`; the density polynomial at
`:1771` is fitted to Earth's T and S range at `sref=35., tref=2.`, and the depth
coordinate uses `du(k)*rhonul/g`. `cpl.f90:3-5` is T21-only AND Earth-radius-only
via `parameter(nxa=64,nya=32)` and `parameter(radea=6.371E6)`. Additionally
`oceanmod.f90:335-344` aborts unless `n_days_per_year == 360`, which this
world's calendar can never satisfy, so LSG is unreachable as the code stands.

**`carbonmod` hardcodes Earth's seconds per year and land fraction.**
Unreachable at `NCARBON = 0`, though `carbonstep` is called every timestep
regardless (`plasim.f90:4088`). `carbonmod.f90:174` uses
`localprecip = (dprc+dprl)*3.154e9`, Earth's seconds per year, two lines from a
`timeweight` at `:114` that correctly uses `m_days_per_year` -- so the two
disagree by a factor of two here. `:194` and `:197` use
`exp(kact*(tsurf - 288.0))`, Earth's global mean as the activation reference.
`:213` divides by `RAD_EARTHSQ*0.29`, where the radius is correctly the model's
but 0.29 is Earth's land fraction against this planet's 50 per cent, and `:27`
hardcodes `RAD_EARTH=6371220.0` beside a correct `plarad`. `tune1 = 5.41` and
`tune2 = 2.20` are commented as "tuning adjustment to make global average
match", so fitted to Earth by construction.

**`icemod_template.f90` is a stale pre-fork copy.** `icemod_template.f90:19`
still carries `parameter(TFREEZE=271.25)`, the `parameter` form from before the
fork made it a namelist key. It has no `nseaice`, no threadprivate directives,
and is not in `CMakeLists.txt`. It is the file a grep for `TFREEZE` hits first
and it gives the pre-CLIM-17 answer. `plasim/bld/make_plasim`, the old unused
makefile, likewise still names `OCEAN=lsgmod` against `lsg/src`.

**`newsnow.f90` and `buildice.f90` are orphaned and hardcode a T21 grid.**
Nothing in `vendor/exoplasim` or `exoplasim/` references either; they are never
compiled and never run. Both hardcode `NLAT=32, NLON=64`, which happens to match
the current rung and would silently mis-read any other. `buildice.f90:21` adds
400 m water equivalent to every land cell unconditionally. They are the other
side of the external mass-balance handshake `glacierstop` writes at
`glaciermod.f90:611-620`, which the model produces on every run and nothing
reads.

---

## Checked and clean

Recorded so they are not re-derived. Several are places where the fork has
already done this work correctly, and knowing which is what makes the findings
above legible as exceptions rather than as the rule.

**Every column amount divides by the model's own gravity.** Ozone
(`radmod.f90:2519`), water vapour (`:2531`), shortwave CO2 (`:2545`), Rayleigh
(`:2553`, `:2761`), longwave H2O/CO2/O3 (`:3249`), CH4/N2O air mass (`:3262`),
cloud liquid water path (`:2465`), and the dust layer thickness. Rayleigh
additionally carries an explicit `9.80665/ga`. Nothing computes a column with
Earth's gravity, and I found no `1/g` omitted anywhere in `rainmod` or
`fluxmod`.

**No Earth boundary field survives into a run.** Codes 169, 210, 211, 709 and
903 all print `Init ... internally` in the diag. `configure()`'s `rm *.sra` does
its job, and the run directory holds exactly the seven project-generated files.
Findings 1 and 19 are what the model does INSTEAD, not a surviving file.

**Every derived scaling is planet-general.** `ww = TWOPI/sidereal_day`,
`acpd = gascon/akap`, `adv`, `cv = plarad*ww`, `ct = cv^2/gascon`,
`rdbrv = gascon/RV`, `pnu21` (`plasim.f90:1514-1521`), and `plavor = EZ`, which
is the non-dimensional planetary vorticity in units of `ww` and correct for any
prograde rotator. `akap = 0.286` with `gascon = 287.017` gives
`acpd = 1003.6 J/kg/K`, which is the right answer for this N2/O2/Ar mixture
rather than an inherited error.

**The planetary constants reach the model correctly.** `PLARAD = 7645464.0`,
`GA = 12.81`, `GASCON = 287.017`, `ECCEN = 0.02`, `OBLIQ = 32.0`,
`GSOL0 = 1286.145`, `ROTSPD = 0.8`, `SIDEREAL_YEAR = 15794043.12`,
`NFIXORB = 1`, `NGENKEPLERIAN = 1`, `PSURF = 100000`, all verified in
`run_2c42c68fe9ca`'s staged namelists rather than in documentation. `meananom0`
presets Earth's Jan-1 value in `p_earth.f90:35` and all three drivers pass
`meananomaly0=0.0`, so it never fires.

**Charnock is gravity-correct.** `seamod.f90:262-263` computes
`z0 = charnock*|tau|*gascon*T/(ga*dp)`, exactly `alpha*u_*^2/g` with `ga`
explicit, and `charnock = 0.018` is dimensionless with no g in it. The floor
`dz0sea = 1.5e-5` is a viscous smooth-flow scale, not a gravity quantity.

**The postprocessor is planet-aware everywhere but one constant.** `radius`,
`gravity` and `gascon` are passed explicitly on both branches
(`__init__.py:1395-1397`, `:1421-1423`). Vertical velocity (`pyburn.py:1869`),
geopotential height (`:2079`), streamfunction (`:1953`) and every vector
transform use them. Finding 17 is the sole exception.

**No Earth calendar reaches any artifact.** The netCDF time axis is
`["time","timestep_of_year","timesteps"]`, raw model timesteps, with no
`calendar` attribute, no epoch and no month names. Seasons are defined by solar
longitude (`analyze_climatology.py:31-36`, `:101-103`) with unequal bin weights
handled by `lib/climatology.py`. Finding 19 is confined to the raw `.srv`
headers and the model log.

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
raised it as a finding and it took source to refute. `icemod.f90:354` reads
`xdt = solar_day / real(ntspd)` and `oceanmod.f90:265` the same for `dtmix`,
which against the module-level `ntspd = 40` and `solar_day = 108744.83` from
`plasim.f90:1546` and `:1498` would give 2718.62 s against the atmosphere's
`deltsec = 2700`. It does not, because both names are rebound on entry:
`iceini`'s dummy arguments are `ktspd` and `psolday` (`icemod.f90:265`), it
assigns `ntspd = ktspd` and `solar_day = psolday` at `:297-298`, and
`seamod.f90:107-109` passes `mtspd` and `day_24hr`. So `xdt = 86400/32 = 2700
s`, identical to `deltsec = day_24hr/mtspd` at `plasim.f90:671`. Two different
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
k25v spectrum at runtime** (`radmod.f90:653-804`) and reach `landmod` through
the namelist, confirmed in the diag echo. The 263.16 K snow aging ramp
(`landmod.f90:395`) is `tmelt - 10` and carries no planet dependence. Findings 7
and 13 concern what is applied ON TOP of correct endmembers.

**The compiled convection scheme has no adjustment timescale to be wrong.** Kuo
is a moisture-convergence closure: `zi` is accumulated over `deltsec2` at
`rainmod.f90:680` and redistributed over the same interval at `:1236-1243`.
There is no CAPE relaxation time, no entrainment rate and no precipitation
efficiency to inherit. There is also no precipitation fall speed anywhere --
`mkrain` moves precipitation through the whole column in one step -- so no Earth
terminal velocity is misapplied.

**The 316 nm shortwave cutoff is reverse-engineered from the Sun and is
nonetheless harmless.** `radmod.f90:449`, `minwavel = 316.036116751`, whose own
comment says it "produces zsolar1=0.517 at Teff=5772 K", zeroes 0.353 per cent
of k25v's flux. Every absorptance is divided by its band's flux share and
multiplied back, so the partition cancels identically; the residue reaches
`zsolar1` itself at 0.57 per cent relative, worth about 0.06 W/m2 and 0.0002 in
surface albedo.

**The Earth air-mass magnification is harmless.** `zm = 35./SQRT(1.+1224.*mu^2)`
(`radmod.f90:2504`) has `sqrt(1224) = 34.99`, so it reduces to 1/mu everywhere
except below mu about 0.03, where the cell receives under 3 per cent of
normal-incidence flux. Under 0.1 W/m2.

**`tfrc`'s timescales cancel their `day_24hr` correctly** and relax over
planetary rotations; the sponge sits on sigma 0.025 and 0.094, so it is
stratospheric rather than boundary-layer drag. `miscmod.f90:76` converts days to
seconds using the planet's own rotation and is the one that is right, though
inconsistent in convention with `plasim.f90`.

**Water properties are not Earth properties.** The Magnus coefficients, `RV`,
`ACPV`, `TMELT`, `sicecap`, `sicediff`, `ALS`, `ALV` and
`L_TIMES_RHOH2O` are properties of the condensable. Only their USE is at issue,
in finding 24.

**The shallow-convection sigma diffusivity converts correctly.**
`rainmod.f90:1543-1544` carries `ga**2/gascon**2`, the right Jacobian for m2/s
into sigma2/s, so only the magnitude of `rkshallow = 10.` is unjustified rather
than the conversion.

**`t0 = 250.0 K`** (`plasimmod.f90:882`) is the semi-implicit reference
temperature, a reference rather than a constraint, and the `nconvtime`
gravity-wave check derives its limit from it correctly with `gascon` and `akap`.

**`utilities.f90` and `restartmod.f90` contain no planetary constants.**
`makeareas` builds fractional spherical cell areas from latitudes only.

## Two things that are neither findings nor clean

**The dynamical non-dimensional timestep is exact by coincidence.**
`plasim.f90:650` and `:673` set `delt = TWOPI/ntspd`, where the correct value is
`deltsec*ww`. These agree only when `ntspd == sidereal_day/deltsec`. At 30 hours
they agree exactly at every rung on the ladder -- 108000/2700 = 40, /1800 = 60,
/1350 = 80 -- so the dynamics and the physics currently run on the same
timestep. Move `rotation_hours` to 26 and `delt` becomes `TWOPI/34` against a
correct `TWOPI/34.667`, and the dynamical core would integrate 2 per cent faster
than the physics, silently. The same lattice is what keeps `solang`'s hour angle
advancing exactly 2*pi per model day with no discontinuity. This wants an
assertion in `check_consistency.py` -- `86400*rotation_hours/24 /
(timestep_minutes*60)` must be an integer -- rather than a source change.

**`N_DAYS_PER_YEAR` survives by ordering.** `__init__.py:2664-2666` computes
`max(int(360.0/rotationperiod/12+0.5),1)*12`, which is 288 for this world:
Earth's 360-day year expressed in Vesper days, not Vesper's year. It is harmless
only because `configure_otherargs` in `run_exoplasim.py:154-157` writes 146
through `otherargs`, which `configure()` applies last. All three drivers pass
that dict, so nobody hits it today. `NSTPW = int(7200//timestep)` at
`__init__.py:2653` is likewise 5 Earth days, giving 160 steps, which is 4 Vesper
days and does not divide 5850; `run_stellar_cycle.py:212-217` works around it
explicitly and the other two drivers do not. The downstream symptom is already
documented and patched in `exoplasim/notes/first-output-bin.md`, so this records
the Earth-day origin rather than a new consequence.

## What this audit did not cover

Stated so the coverage is not overread. The `puma`, `cat` and `sam` trees were
not audited, nor the GUI path, the FFT modules, the MPI decomposition, or the
`postprocessor` and `pRT` subtrees beyond the constants they expose. Within the
compiled set, `rainmod`'s cloud-cover diagnostic feeding `dcc` and the leapfrog
`deltsec` against `deltsec2` question were looked at and left; the second is a
water-budget question for `water-and-energy-closure.md` rather than an
Earth-centrism one. Physics that is simply ABSENT was out of scope and is the
subject of `absent-and-inherited-physics.md`.

This audit read source and staged artifacts. It ran no model. Every magnitude
quoted is an estimate from the code and the configuration rather than a measured
perturbation, and the ones most worth measuring are finding 1's sweep behaviour,
finding 2's factor of 1.25, finding 3's factor of 1.44 against the arms already
run, and finding 6's emissivity bracket.

## Tasks

Tracked in the `bd` issue tracker, not restated here. Rows carrying the label
`audit:model-earth-centrism` were opened by this audit; the rest were already
open when it ran, from the sweeps of the same day, and this audit corroborates
them rather than duplicating them. Where two rows cover one finding, both are
named because they were filed independently.

| finding | id |
| --- | --- |
| 0. continuations revert hyperdiffusion and filter power | `world-8bs` |
| 1. the cold start builds a CCM3 Earth Arctic ice field | `world-6fh` |
| 2. the damping timescales are applied 25 per cent too weak | `world-rt1` |
| 3. the ocean's diffusion runs on Earth's radius | `world-st4`, `world-mll` |
| 4. a 5 m snow cap neuters the glacier module | `world-cwc` |
| 5. cloud liquid water uses an Earth scale height | `world-ofn` |
| 6. land longwave emissivity is exactly 1.0 | `world-qvu` |
| 7. the snow-over-forest albedo is one factor for both bands | `world-nfh` |
| 8. the critical relative humidity is an Earth grid-box tuning | `world-khn` |
| 9. precipitation re-evaporation is Earth-calibrated and timestep-dependent | `world-j6v`, `world-khn` |
| 10. salinity is four constants and `CPS` is fresh water's | `world-9hb` |
| 11. the ozone profile sits at a fixed geometric height | `world-ayx` |
| 12. the sea-ice albedo ramp keeps an Earth slope | `world-cj4` |
| 13. the free-convection ocean coefficient embeds Earth's gravity | `world-e2k` |
| 14. one Earth-average soil thermal pair | `world-sy9` |
| 15. the archived albedo 175 is the Earth-Sun broadband value | `world-e5p` |
| 16. mean sea-level pressure uses Earth's lapse rate | `world-ld1` |
| 17. Earth's lapse rate sets the cold-start profile | `world-wmw` |
| 18. the model's calendar is a 360-day Earth year | `world-x1k` |
| 19. glacier persistence is one run segment | `world-qpe` |
| 20. an undeclared Rayleigh sponge of 20 and 100 days | `world-aee` |
| 21. the mixing length is fixed metres in a sigma scheme | `world-e2k` |
| 22. snow masks the surface at 0.01 m | `world-khn` |
| 23. convective cloud cover is an Earth mm/day regression | `world-khn` |
| 24. no saturation-over-ice branch | `world-ako` |
| 25. snow and glacier densities are Earth compaction values | `world-1pl` |
| 26. the lead-closing scale and the thickness clamp | `world-12c` |
| 27. runoff velocity constants | `world-529` |
| 28. the run log identifies the planet as Earth | `world-1o4` |
| mechanism I, the build compiles `p_earth` | `world-cwu`, `world-58v` |
| dormant, dust `APART` and `RHOP` are never written | `world-906` |
| dormant, SIMBA and its two traps | `world-9hv` |
| dormant, `nfluko` relaxes toward the constructed field | `world-4ba` |
| dormant, the radiation tuning is pinned at T21 | `world-ys9` |
| dormant, the rest, with what each waits on | `world-9d1` |
| dormant, deleting what no build compiles | `world-cmz`, `world-6ak` |
| dormant, `newsnow.f90` and `buildice.f90` | `clim-53`, `world-zsa` |

`world-cwu` records `alr` as inert. Finding 17 corrects that, and the correction
is appended to the row.
