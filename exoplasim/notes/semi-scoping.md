# SEMI as the surface mass balance candidate: what it needs, and what it brings

Worldbuilding frame: this is a scoping read of external source against the
Vesper project's climate model. Every field named is a model array or a model
input; nothing here is about the real world.

*Read 2026-08-24 from `references/climber-x/src/smb/`, which is READABLE
REFERENCE and not a vendored component. Nothing was adopted in this pass: no
`vendor/` addition, no `config/pipeline.yaml` step.*

SEMI is CLIMBER-X's surface energy and mass balance model. `clim-63` asks
whether it is the candidate and what it would cost; `dust-14` waits on the
answer because the user's decision was to solve dust-on-snow once, correctly,
inside whatever brings a snow grain size, rather than twice.

## 1. The grain size arrives with SEMI, and it arrives coupled to the dust term

This is the question the decision to wait turns on, and the answer is yes.

`smb_surface_par.f90` is one module and it holds all three pieces. `surface_albedo`
calls `snow_grain_size` when `snow_par%lsnow_aging` is set, calls `dust_in_snow`
when `snow_par%lsnow_dust` is set, and then hands BOTH results to the same snow
albedo routine. The two are arguments to one function, not two independent
corrections:

    call snow_albedo_dang(z_sur_std, snow_grain, dust_con, coszm, &
                          alb_snow_vis_dir, alb_snow_nir_dir, alb_snow_vis_dif, alb_snow_nir_dif)

Inside `snow_albedo_dang`, the coupling is explicit rather than incidental. The
grain radius sets the clean-snow albedo through `rn = log10(r/r0)` with `r0` at
100 um, and the same radius then scales the darkening: the black-carbon-equivalent
loading is `H = c/c0*(r/r0)**0.73`, so a dust concentration darkens coarse snow
more than fine snow, by a stated power. That exponent IS the reason a dust
coefficient cannot be transplanted onto a grain-free ramp. `snow_albedo_ww`, the
alternative under `isnow_albedo == 1`, reaches the same structure by a different
route: it interpolates a tabulated Warren and Wiscombe dust curve in
`log10(dust_con)` and then weights it by a snow-age factor derived from the same
`snow_grain`, with separate coefficients for new and aged snow.

So the premise the decision rested on holds. Adopting SEMI settles dust-14's
blocker, because the grain size, the dust concentration, and the calibration that
relates them are one object.

**The calibration that comes with them**, from `references/climber-x/nml/smb_par.nml`
and the `snow_par_type` defaults in `smb_params.f90`:

- `isnow_albedo` selects Dang et al. 2015 over the CLIMBER-2 Warren and
  Wiscombe 1980 scheme, and both are shipped.
- The grain size runs between a fresh value and an old value, both declared, and
  the ageing law between them is a CLIMBER-2 parameterisation tuned to MARv3.6
  using the CROCUS snow model over Greenland. That provenance is stated in a
  comment in `snow_grain_size` and it is an Earth calibration of a snow process,
  not of a stellar spectrum.
- `dust_con_scale` exists precisely to rescale the dust concentration for a
  different imaginary refractive index, which is the knob a different dust
  mineralogy would turn. It defaults to unity, so the shipped state is Dang's
  own calibration untouched.
- `w_snow_dust` is the melt concentration length: the snow water equivalent whose
  melt doubles the dust concentration, capped at a five-fold increase, on the
  reasoning that meltwater scavenges only 10 to 30 percent of the dust.

**The dust term needs no prognostic reservoir.** `dust_in_snow` computes
`dust_con = dust_dep/max(1e-7, snow)`, a concentration in FALLING snow, times the
melt factor above. The only state involved is `w_snow_max`, the seasonal maximum
snow water equivalent, one scalar per cell that `snow_update` raises whenever the
pack is accumulating. `notes/audits/absent-and-inherited-physics.md` finding 2
concluded that closing this coupling costs "a prognostic snow-dust reservoir, its
restart record, an albedo function that reads it". The peer's answer is cheaper
than that: a flux-ratio diagnostic, one extra scalar of state, and the albedo
function.

**And this project already holds the grain-resolved reflectance data**, which
turns out to matter more than it first looks.
`vendor/exoplasim/exoplasim/plasim/src/specblock.f90` ships `fsnowalb`,
`msnowalb` and `csnowalb`: fine, medium and coarse granular snow reflectance over
965 wavelengths, from the same JHU becknic library every other endmember here
came from. No code path selects among them -- `radmod.f90` uses the `iceblend`
family instead, and the one line that would have combined the three is commented
out. `references/exocam/tools/spectral_albedos/snow100um.txt` is a second anchor.

The `iceblend` family is not a grain axis, and it is worth saying so precisely,
because it looks like one. `plasim/src/specs/combinedspec-snow_ice.ipynb` is the
notebook that built those arrays, and it is in the tree. Each blend is a weighted
mixture of five component spectra -- clear ice, frost, and coarse, fine and
medium granular snow -- with the snow components at fixed proportions and the
CLEAR ICE fraction solved by bisection against a declared broadband albedo target
under a 5772 K blackbody. So `iceblendmin` and `iceblendmax`, the endpoints
`landmod`'s snow albedo ramp interpolates between on surface temperature, differ
in how much clear ice is mixed in, tuned to hit a solar broadband number. They do
not differ in grain size, and the temperature ramp between them is not an ageing
law in any physical sense.

That is why a dust coefficient grafted onto that ramp would be doubly unanchored:
not only is there no grain radius for the darkening to scale with, the two
endpoints it interpolates between are themselves mixture fits to an Earth target.

It also points at the constructive route, and section 6 returns to it. The three
grain-resolved spectra plus `analysis/ice_albedo.py`, which already reproduces
`radmod`'s star-weighting arithmetic to machine precision, are between them
enough to derive a two-band snow albedo as a function of grain radius FOR THIS
STAR, on a three-point axis, from data already in the tree. That is the half a
dust coefficient needs in order to mean what its source measured, and it does not
require SEMI.

## 2. The four-component structure is real, but it is not direct and diffuse

`clim-63` and `notes/external-model-survey.md` section 3f both describe SEMI's
albedo and shortwave inputs as "visible and near-infrared by direct and diffuse".
The count and the spectral axis are right; the second axis is not. The suffixes
`_dir` and `_dif` in this module mean CLEAR SKY and CLOUDY SKY, and
`smb_def.f90`'s own comments say so:

    swd_sur_vis_dir   !! downward shortwave visible radiation at surface, clear sky [W/m2]
    swd_sur_vis_dif   !! downward shortwave visible radiation at surface, cloudy [W/m2]

`rad_downscaling` then recombines them by cloud fraction, `(1-cld)*clear +
cld*cloudy`, rather than by a beam-splitting ratio, and `coupler.f90` wires
`dswd_dalb_vis_dir` from the atmosphere's `dswd_dalb_vu_cs` -- `cs` for clear
sky. A direct-beam versus diffuse-beam reading would send this project looking
for a decomposition ExoPlaSim does not make and does not need to make.

The zenith angle enters separately, through `coszm`, in the snow albedo itself.

**The spectral axis IS this project's.** `constants.f90` sets `frac_vu = 0.45`,
the fraction of the solar spectrum in the visible and ultraviolet, and
`rad_downscaling` uses it to weight the two bands into a broadband flux. That is
the same quantity `lib/stellar.py` computes at the 0.75 um split and
`world_state.json` records, so the substitution is one constant and it is one
this project has already derived for this star. The mapping in section 3f is
verified.

## 3. Field by field: what ExoPlaSim emits, what is derivable, and what is absent

SEMI's column takes 32 inputs per cell per step. Grouped by what this project
would have to do:

**Already emitted.** Total precipitation, total cloud cover, surface wind speed,
surface soil temperature, and the two-band surface albedo. The albedo is the
useful one: `dsalb(1,:)` and `dsalb(2,:)` are archived as codes 174 and 184, with
the broadband as 175, and after the wave-1 surface albedo work `landmod`'s snow
albedo endpoints are themselves per band, `albsmin1/albsmax1` and
`albsmin2/albsmax2`, driven from `dsnowalbmn`/`dsnowalbmx` which `radmod`
computes by integrating a measured reflectance blend against this star's
spectrum. SEMI wants four albedo components and ExoPlaSim has two: its surface
albedo carries no cloud dependence, so the clear and cloudy arms are the same
number, and the mapping is two fields used twice rather than two fields missing.

**Derivable from what a run already writes.** Free-atmosphere temperature and the
near-surface lapse rate, from the temperature profile and `lib/lapse.py`, which
already owns the height of the lowest model level. Relative humidity, from the
specific humidity and temperature. The daily-mean cosine of the zenith angle,
from the orbit. Downward longwave at the surface and downward shortwave at the
surface exist as `dftd(:,NLEP)` and `dfd(:,NLEP)`, codes 407 and 405, but only
inside the level-resolved output rather than as surface diagnostics -- codes 176
and 177 are the NET surface fluxes, not the downward ones.

**Derivable only from a run configured to produce it.** The standard deviation of
daily 2 m temperature, `tstd`, which `smb_ebal` uses for the sub-diurnal melt
correction, needs daily sampling within the averaging period. The diurnal minimum
of top-of-atmosphere downward shortwave, `swd_toa_min`, is the same class. The
700 hPa winds are available through the postprocessor's pressure-level
interpolation, but 700 hPa is itself an Earth choice -- it is the steering level
for Greenland orographic precipitation -- and on a planet with a different scale
height it is a declared choice rather than an inherited one.

**Not emitted, and this is the model change.** The band-resolved downward
shortwave at the surface, in both the clear-sky and the all-cloud limit: four
fields. `radmod`'s `swr` computes `zfd1` and `zfd2` per band at every interface
and then immediately sums them into `dfd`, so the band split exists inside the
routine and is discarded. The clear-sky arm has a precedent already in the tree:
`ndiagcf` re-calls `swr` with the cloud fraction zeroed and stores the result in
`dclforc`. The all-cloud arm has no precedent -- CLIMBER-X's atmosphere computes
`solar_sur_c` alongside `solar_sur_s`, an overcast limit rather than the actual
sky, and ExoPlaSim's normal call gives the cloud-weighted sky only.

**Not derivable in closed form.** The six derivatives: `dswd_dalb` for each of
the four components, and `dswd_dz_nir` for the two near-infrared ones. In
CLIMBER-X these are analytic derivatives of that model's own single-layer
transmission-factor product, written out by hand in `swr.f90` around the
`l_dswd_dalb` guard. They are inseparable from that scheme. ExoPlaSim's shortwave
is a multi-layer adding cascade with the surface albedo entering as the bottom
boundary reflectivity, so the derivative would have to be re-derived through the
cascade or taken by finite difference -- two extra `swr` calls per radiation
step, on the `ndiagcf` pattern. Either is real work and neither is a port.

**Deliberately dropped.** `t2m_bias`, `prc_bias` and everything
`smb_bias_corr.f90` does are corrections against Earth observations, of which
this world has none. `dTvar` is an artificial interannual variability generator.
`l_regional_climate_forcing` reads a regional climate model's output instead of
downscaling, which is not a mode that exists here.

## 4. The sub-grid side, and one input this project already has

`clim-63` and section 3f list "sub-grid elevation, an elevation standard
deviation, slopes and elevation classes". Three of those four are right and the
fourth is not a SEMI structure at all.

**Elevation classes do not exist in SEMI.** The only occurrence in the module is
in `smb_out.f90`, binning ice AREA by elevation for output. SEMI does not tile a
cell into elevation bands. It runs the whole balance on a separate, finer grid --
`ice_grids.nml` defines polar stereographic domains from 80 km down, centred at
70 degrees north and 45 degrees west, which is Greenland. So the sub-grid
structure SEMI needs is a finer GRID plus per-cell statistics on it, not a
sub-grid tiling.

**The elevation standard deviation exists here today.** The Orogen export already
carries `orog_std` per grid cell, with `orog_mean`, `orog_min`, `orog_max`,
`orog_count`, `orog_anisotropy` and `orog_angle` beside it, in `planet.nc` and in
`grid/`, on every emitted grid. The manifest's `subgridOrography` block documents
them and states what an empty cell means. That is `z_sur_std` directly: SEMI uses
it for the snow cover fraction over rough topography, for the ice fraction in
`ice.f90`, and for a roughness albedo reduction in `snow_albedo_dang`. It is worth
saying plainly because `grid-2` reads as though no sub-grid elevation statistic is
persisted anywhere, and this one is. What `grid-2` is actually about -- the
DISTRIBUTION of mesh elevations under a cell, which mean, standard deviation,
minimum and maximum do not determine -- is a different and larger object, and
SEMI does not need it.

CLIMBER-X derives its own `z_sur_std` in `geo/hires_to_lowres.f90` by the law of
total variance, combining a read-in sub-grid variance with the variance of the
high-resolution cell means. The export's `orog_std` is the second of those two
terms only, computed over mesh region centres. `notes/audits/orogen-resolution.md`
already establishes that the generator designs nothing below about 20 km, so the
first term is not available here at any region count and is `gw-6`'s.

**The slopes would have to be built.** `topo_grad_map` needs `dz_dx_sur` and
`dz_dy_sur` on the fine grid, for the wind-slope precipitation factor. The export
gives `orog_anisotropy` and `orog_angle`, the eigenstructure of the slope
covariance tensor within a coarse cell, which is a different object: it describes
whether the sub-grid roughness is ridged and along which axis, not the resolved
gradient of the fine field. Slopes on whatever fine grid is chosen would come
from the mesh, and the gradient formula in `topo_grad_map` hardcodes `r_earth`
and must not be ported as written.

**The fine grid itself does not exist.** This is the real gap, and it is the
same one `phys-13` and `grid-2` describe from two other directions. SEMI needs a
grid on which land is resolved finely enough that the balance is not evaluated on
a cell mean, an ice fraction and an ice albedo on it, and a mapping to and from
the climate grid. `phys-13` establishes that the model's own orography evaluates
its high ground far too warm for ice to survive, and that the sub-grid peak excess
does not go away at any truncation the model can afford. A mass balance run on
the coarse orography would return no ice, confidently and for the wrong reason,
which is `clim-53`'s stated ordering argument and it applies to SEMI unchanged.

## 5. Does SEMI survive the test that refused the positive-degree-day scheme?

`clim-63` refuses MITgcmIS's PDD scheme for driving ablation from 2 m air
temperature alone, blind to absorbed shortwave over ice, which is the term a
K2.5V host changes most. SEMI is checked against the same test rather than
assumed to pass it.

**On the shortwave, it passes cleanly.** `smb_ebal` linearises a full surface
energy balance and takes the melt flux as the residual at the melting point:
`flx_melt = sh + lh + g + lw + swnet` evaluated with the skin temperature pinned
to freezing, floored at zero. `num_sw = swnet` enters the skin temperature solve
directly. `swnet` itself is assembled in `rad_downscaling` band by band against
the band-resolved surface albedo. There is also a sub-diurnal correction that
uses `swnet_min` to estimate the fraction of the step spent above freezing, so
the scheme resolves melt that a step-mean temperature would miss. Absorbed
shortwave over ice is the central term, not an absent one. It also ships a
surface energy conservation check, `energy_cons_surf1` under `check_energy`,
which is the kind of identity CLAUDE.md's testing convention asks for and which a
port could keep as its acceptance test.

**On the snow albedo, it is a step backwards in spectral fidelity, and that is
the one thing about SEMI worth arguing over.** Dang et al. 2015's coefficients --
0.9856, -0.0202, -0.0125 for the diffuse visible clean-snow albedo, and the
parallel sets for the other three components -- are fits to band-INTEGRATED
albedo under a solar spectrum, with a visible/near-infrared boundary that is that
scheme's own and not 0.75 um. Under a K2.5V host, more flux sits in the near
infrared and the in-band spectral distribution shifts, so the band-integrated fit
no longer describes the band. Dang's polynomials cannot be re-derived here,
because the scheme ships the fits and not the spectra they came from.

ExoPlaSim's endmembers are in the opposite position, and `phys-14`'s overturned
verdict is the evidence: `radmod` already integrates the shipped reflectance
blends against the configured stellar spectrum, splits at 0.75 um, and assigns
this star's numbers rather than solar defaults. `analysis/ice_albedo.py`
reproduces that arithmetic as the registered `ice_albedo` step. So on the axis
that matters for a non-solar host, this model is currently AHEAD of the scheme it
would be adopting, and behind it on grain size and impurities. The two schemes
are strong in different places.

So SEMI survives the refusal that killed the PDD scheme. An energy balance with a
spectrally imperfect albedo is strictly better than a temperature index with no
shortwave at all, and physics is not a knob. But the right adoption is not a
straight port of the albedo half, and section 6 says what it is instead.

**The other inherited constants**, all PHYS-class, in one list so a port does not
discover them one at a time: `frac_vu = 0.45` (section 2, this project has the
replacement); `p0 = 1010 hPa` and `h_atm = 8600 m` in `smb_params`, both carrying
`fixme` comments upstream, used by `topo_factors` to make a surface pressure --
the scale height is gravity-dependent and must come from `config/planet.yaml`;
`gamma = 5 K/km`, a declared lapse rate with an Earth ice-sheet provenance, where
`lib/lapse.py` owns this project's; `r_earth` in `topo_grad_map`; `g` in the
Richardson number in `resistance` and in `snow.f90`; the roughness lengths
`z0m_ice` and `z0m_snow`; the emissivities; `dP_dT`, the Clausius-Clapeyron
precipitation factor; and `wind_ele_fac`, whose own comment says to check it
against CORDEX data. Snow density is a constant here, so unlike PALADYN's snow
there is no compaction law and section 3b's gravity argument does not apply to
this module.

## 6. The adoption shape: a port, not a subtree, and not yet

`docs/src/reference/vendored-upstreams.md` sets the rule. A subtree is a
MAINTAINED FORK, taken when the work is fork-shaped -- when this project will
edit the external code and wants provenance as a commit here. External source
this project reads and will not edit lives under `references/` instead.
`vendor/cgenie` is the boundary case and it is instructive: it is a subtree
because two open rows are fork-shaped by construction, and it is still labelled a
CANDIDATE that nothing reads.

SEMI does not meet that bar, and the reason is that almost nothing of the
surrounding tree transfers.

The module is 8,378 lines. The COLUMN -- `semi.f90`, `smb_surface_par.f90`,
`smb_ebal.f90`, `smb_temp.f90`, `snow.f90`, `smb_grid.f90`, `downscaling.f90` --
is about 1,900 of them, and its dependencies outside `smb/` are small and
ordinary: a precision kind, a physical constants module, a namelist reader, a
tridiagonal solver. That part is a clean, self-contained energy balance on a
snow-over-ice column and it is genuinely portable.

The other 6,400 lines are the driver, the grid machinery, the Earth bias
corrections, the polar stereographic Greenland and Antarctic domains, the netCDF
output, and CLIMBER-X's own timer and coupler contracts. None of it transfers.
`smb_bias_corr.f90` has no meaning on a world with no observations, and
`smb_simple.f90` and `smb_pdd.f90` are the alternatives `clim-63` already
refused.

So the shape is a PORT of the column into this project's own component
directory, with the parameter block re-derived rather than copied, against the
existing components -- not a subtree, and not a reimplementation from the paper
either, because the code is the specification and the licences combine.
GPL-3-or-later over ExoPlaSim's GPL-2-or-later is permitted with the combined
work under GPL-3, which section 3d of the survey establishes.

**Two things must exist before a port is worth starting**, and both are other
rows:

1. A fine grid with an ice mask on it. `phys-13`. Until then the balance runs on
   an orography that evaluates its own high ground far too warm, and returns no
   ice for the wrong reason.
2. The four band-resolved downward shortwave surface diagnostics, in the
   clear-sky and all-cloud limits, and a decision on how the six albedo and
   elevation derivatives are obtained. That is a change to `radmod`, `outmod`,
   `pyburn` and the diagnostic block, and it is the half of `clim-63` that is
   ExoPlaSim work rather than CLIMBER-X reading.

**And one piece is worth doing FIRST, independently of both.** Section 5 leaves
the albedo half of SEMI as the part that should not be ported as written. The
replacement is derivable from data already in the tree: take the three
grain-resolved snow spectra out of `specblock.f90`, weight each against this
star's spectrum through `analysis/ice_albedo.py`'s existing arithmetic, and get a
two-band snow albedo on a three-point grain-radius axis for this star. Fit Dang's
functional form to THOSE points rather than importing his coefficients, and the
dust term's grain scaling then has a local anchor.

That is a small, self-contained piece of work with a check that can fail -- the
solar-weighted arm must reproduce the published values within the tolerance
`ice_albedo.py` already establishes for the existing blends -- and it does not
need SEMI, a fine grid, or `phys-13`. It is also the shortest honest route to
closing `dust-14`, because it supplies exactly the missing half: a grain size for
a published dust coefficient to scale against. Whether the grain size then
EVOLVES, and by what law, is the part that still arrives with a mass balance.

Item 2 has a cheaper first step that is worth naming: SEMI's derivatives exist to
let the balance re-solve the surface albedo on the fine grid without re-running
the radiation. A first port could hold the albedo fixed at the climate model's
value, drop all six derivatives, and lose only the albedo feedback WITHIN the
downscaling step. That is a defensible reduced form with a stated invalidity --
it cannot represent a snow line moving inside a coarse cell during the melt
season -- and it removes the hardest of the input requirements from the critical
path.

## 7. What a mass balance is allowed to claim, given no time axis

`docs/src/reference/no-time-axis.md` binds anything that asks for a duration.
SEMI itself is clean: it produces a rate in kg/m2/s, and every duration in it
comes from the climate -- the model year, the SMB timestep, the climate update
interval -- not from the terrain. Orogen is not asked for anything.

Where it bites is the step after. A surface mass balance field is not an ice
sheet. Turning one into ice geometry needs an ice dynamics model integrated over
a duration, and there is no duration to integrate over; a declared window would be
exactly the failure the document names, a number that is calibrated to an outcome
being reasoned about as though it were measured.

The document's own method gives the way through, and it is the standard one: ask
what a randomly chosen moment looks like, which for a mass balance means solving
for the EQUILIBRIUM rather than a history. The honest deliverables are the
annual balance field itself and the surface at which it integrates to zero. Both
are diagnostics of the climate, and neither needs a duration. Anything that asks
how THICK the ice is, or how long it took, does, and `clim-53` is where the
acceleration question and its declared invalidity already live.

That also fixes what `dust-14` gets. The dust term changes the melt energy and
therefore the equilibrium line, which is a statement the world can carry. It does
not by itself say how much ice there is.

## 8. The criteria that decide the adoption, and the question sections 1 to 7 do not ask

Declared 2026-09-05, before the cost below was worked out. Sections 1 to 7 settle
that SEMI is the candidate on its STRUCTURE. What CLIM-63 still owes is a price
and a criterion, and one question those sections do not ask.

**Ask first whether the thing should exist.** Deleting is a fix, so the port has
to beat two alternatives and not one. The refused alternative is the
positive-degree-day scheme, and section 5 disposes of it. The alternative that
has never been named is **the model's own land surface**: `landmod` already runs
a surface energy balance with snow, carries a snow water equivalent, and melts
it. If what SEMI adds over that is only the downscaling, the grain size and the
dust term, then the port is those three pieces and not a mass balance scheme, and
the price is smaller than section 6's 1,900 lines. That question is answered
below before any cost is quoted.

**The criteria.**

- **SEMI IS ESTABLISHED as the candidate** if it survives the test that refused
  the PDD scheme (section 5: it does), if every input it needs is either
  available now or held by a named row, and if no cheaper object in the tree
  already computes what it computes.
- **THE PORT IS BUYABLE** if its recurring run cost is small beside a
  commissioning. A surface mass balance runs offline on a climatology rather than
  inside the climate model, so the bar is that one evaluation over the fine grid
  costs less than one ORBIT of the climate model at the configured rung -- 13.45 s
  at T21, 75.74 s at T42, 491.96 s at T85, from
  `notes/audits/resolution-ladder-wall-clock.md`. Above that the mass balance
  becomes a term in loop A's cadence rather than a diagnostic of it.
- **AND THE PORT IS REFUSED, whatever it costs**, if what it would deliver needs
  a duration. Section 7 already binds this: the honest deliverables are the annual
  balance field and the surface at which it integrates to zero.

**Rule 7.** Nothing is commissioned, so nothing is charged for invalidation. The
cost of adopting SEMI is what building it costs plus what the next cycle runs.

## 9. What the model already has, read 2026-09-05, and it is not a mass balance

The alternative section 8 says has never been named is `landmod`'s own snow. It
was read, and it does not do the job. What exists is a one-layer bucket snow
store with an energy-balance melt term; what a mass balance needs beyond that is
absent in five separate places.

**The melt term is the good half and it is the right shape.** `tands` at
`landmod.f90:1763-1781` clamps the surface temperature at the melting point and
takes the melt flux as the residual of the surface energy balance against that
clamp, with the atmospheric forcing `zhfla = dshfl + dlhfl + dflux(:,NLEP)`.
That is structurally SEMI's `smb_ebal` idea, and it is why the positive-degree-day
refusal in section 5 was never about this model.

**What is absent.**

- **Cold content, because there is no snow temperature.** `dsnowt` is declared,
  carried in the restart, and set in four places -- `landmod.f90:686`, `:1704`,
  `:1837`, `:2366`. Nothing anywhere READS its value into a flux or a temperature.
  The model's snow is isothermal with the surface skin, and the source says so at
  `:1834`. A write-only prognostic in the restart is its own defect and is filed
  separately; what it costs HERE is that there is no physical basis for the next
  item.
- **Refreezing, entirely.** Melt goes into `dwater`, the soil-water forcing, and
  from there to the bucket and to runoff. `landphase` freezes and thaws SOIL
  water against the soil temperature; nothing returns water to the pack.
  SEMI's `smb_temp.f90` charges the cold content before the latent heat of
  fusion and splits the result into snow and superimposed ice.
- **Rain on snow.** `dprl` and `dprc` go straight to `dwater` and never reach the
  pack.
- **Ice melt as a term distinct from snow melt**, and therefore an ablation zone.
  There is only `dsmelt`.
- **Any accumulator at all.** Nothing in `outmod.f90`, in `pyburn.py`'s code
  table or in `exoplasim/scripts/` computes accumulation minus ablation, an
  annual net balance, or an equilibrium line. The nearest thing is `asndch`,
  which accumulates the per-step change in `dsnowz` without normalisation, and
  it is a difference of the store rather than a balance.

**And one term is wrong rather than missing.** When `dsnowz > 0` the WHOLE cell's
`devap` is charged to the snowpack, and `fluxmod.f90:913-921` gives it the LIQUID
latent heat wherever the skin is above freezing. SEMI's `update_tskin` diagnoses
a snow-specific sublimation over ice saturation instead.

**The clamp does not bite here and it is worth recording why.** `dsmax` routes
snow above its ceiling into `dwater` and reports it as `dsmelt`, so under a
positive `dsmax` the melt diagnostic carries numerical discard. Every arm in
`exoplasim/analysis/arms/` sets `max_snow_depth_m: -1.0`, so this model
accumulates without a ceiling and `glaciermod` converts persistent deep snow to
glacier on a threshold. That threshold -- `dsnowz` above `glacelim` for
`glacpersist` orbits -- is a persistence rule and not a balance either.

**So the answer to section 8's question is that SEMI is not duplicating anything.**
The overlap is one term of eight, and it is the term the model already gets right.

## 10. The price, and the half of it that does not wait for the fine grid

**The port is larger than section 6 says.** The seven modules that section names
are 2,049 lines, not about 1,900: `semi.f90` 292, `smb_surface_par.f90` 446,
`smb_ebal.f90` 431, `smb_temp.f90` 446, `snow.f90` 91, `smb_grid.f90` 81, and
`downscaling.f90` 262, the last carrying module `downscaling_mod` under a
different file name. Two things section 6 does not count come with them:
`smb_params.f90` at 370 lines, which six of the seven `use` and which pulls
`ncio`, `nml` and `timer` unless it is rewritten against this model's own
namelist machinery; and `tridiag.f90` at 133 lines, which `smb_temp` needs and
which this model has no exposed equivalent of, since `mktsoil` hand-rolls its own
soil solve. **Call it about 2,550 lines touched, of which roughly 1,700 is
physics this model does not have.** `precision`, `control` and `timer` are
substitutions rather than ports, and `constants` must be re-pointed rather than
copied because this model already owns its saturation functions.

**AND THE PORT SPLITS CLEANLY IN TWO, WHICH THE ORDERING ARGUMENT DID NOT
ANTICIPATE.** CLIM-53's argument, which put this row after PHYS-13, is that a
mass balance fed cell-mean orography answers the wrong question. That argument
binds the DOWNSCALED half and not the column half.

- **The column half is `smb_temp.f90`, `snow.f90`, `smb_grid.f90`, `tridiag.f90`
  and the annual accumulators: under 800 lines.** It supplies cold content,
  refreezing, rain on snow, ice melt as its own term, and the annual balance
  accumulators. Every one of those is a defect in the model's own column, on the
  grid it already runs, and none of them needs a fine grid, a band-resolved
  surface shortwave, or a snow albedo rewrite. It does not depend on PHYS-13 and
  it does not depend on WORLD-QK4I.
- **The downscaled half is `downscaling.f90`, `smb_ebal`'s diurnal branch and
  `smb_surface_par.f90`: about 1,750 lines.** This is the half CLIM-53's
  ordering argument is about, and it is the half that needs the fine grid, the
  four band-resolved downward shortwave fields in both sky limits, and the
  fine-grid slopes that `orog_anisotropy` and `orog_angle` are not.

**The sub-diurnal correction belongs to the second half and is worth naming**,
because it is the largest thing in SEMI that this model structurally cannot
reproduce. `smb_ebal.f90:154-244` reconstructs the diurnal skin-temperature
amplitude from `swnet` and `swnet_min` -- the daily minimum net surface
shortwave -- adds a synoptic standard deviation, computes analytically the
fraction of the day above freezing, and takes the melt residual at the mean
temperature over that fraction alone, then charges the extra melt back against
the ground heat flux as overnight refreezing. It exists to let a cell whose
DAILY MEAN skin temperature is below freezing still melt. This model's land step
already runs sub-hourly, so on the coarse grid it resolves what that correction
reconstructs; the correction earns its place only when SEMI is driven from
downscaled daily or monthly forcing on a fine grid, which is the intended use.

**The run cost is not the price.** SEMI is a per-column tridiagonal solve over a
fixed layer count plus closed-form surface algebra, and the criterion fixed in
section 8 was that one evaluation over the grid must cost less than one orbit of
the climate model at the configured rung. A tridiagonal solve of ten to twenty
layers over a few thousand columns is milliseconds against 13.45 seconds an
orbit at T21, so the column half clears that criterion by orders of magnitude and
the question does not need measuring. The downscaled half's cost is set by the
size of a fine grid that does not exist, and cannot be quoted until PHYS-13
declares one.

**And section 7's limit is unchanged and binds both halves.** What either half
may honestly deliver is the annual balance field and the surface at which it
integrates to zero. An ice thickness needs a duration and there is none.

## 11. The verdict

**SEMI IS ESTABLISHED as the surface mass balance candidate.** It survives the
test that refused the positive-degree-day scheme, its inputs are each either
available now or held by a named row, and nothing in this model already computes
what it computes -- the overlap is one term of eight and it is the term this
model already gets right.

**The column half is BUYABLE NOW and is not blocked.** Under 800 lines, no fine
grid, no band-resolved surface shortwave, no albedo rewrite, and it closes five
absent terms in a column this model runs today. CLIM-53's ordering argument does
not reach it.

**The downscaled half stays ordered behind PHYS-13**, on that argument unchanged,
and its price cannot be quoted until a fine grid is declared.

**Rule 7.** Nothing is commissioned, so neither half is charged for output it
would invalidate. The price of each is what building it costs plus what the next
cycle runs.
