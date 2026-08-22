# What other models in this class do, and what is worth taking

*Written 2026-08-22. Claims about this tree are from reading it; claims about
each external model are from its published description or its source, cited
inline. Costs quoted for other models are what their own authors state, on their
own hardware, and are not benchmarks run here.*

Vesper is a fictional planet and this project simulates it. What follows
compares this project's machinery against the Earth system models of
intermediate complexity, the exoplanet general circulation models, and two
recent projects that sit between them, so the comparison is not run again. It
reaches no adoption decision. What it produces is a short list of formulations
worth taking and a shorter list of candidates that are now closed.

## 1. The two families do not overlap, and this stack is in the gap

Every EMIC in current use runs one planet. CLIMBER-X, CLIMBER-2, LOVECLIM, UVic
ESCM, Bern3D, cGENIE, MIROC-lite and DCESS are Earth-configured down to their
constants: none accepts a non-solar host spectrum, a non-Earth radius, or a
rotation rate outside a narrow tuning range. Paleogeography is the most any of
them varies.

Every model that does run other planets is an atmosphere with almost nothing
under it. ExoCAM, ROCKE-3D, LMD-G, the Met Office UM, Isca, Exo-FMS and stock
ExoPlaSim ship no carbon cycle, no soil model, no drainage network and no
dynamic global vegetation model between them.

ExoPlaSim is the only entry with a foot on both sides, and only just: its
glacier and carbon-silicate weathering modules exist, are off here, and carry
Earth constants when on. `notes/audits/dormant-exoplasim-modules.md` records
which.

## 2. Cost, in CPU-hours per 100 simulated years

The unit matters. Simulated years per day hides the core count, and these models
were measured on different machines with different core counts by different
people.

| model | configuration | CPU-h per 100 sim-yr | source |
| --- | --- | ---: | --- |
| CLIMBER-X | 5 x 5 deg, full climate component | ~3.8 | 10,000 sim-yr/day on 16 CPUs, Willeit et al. (2022) |
| LOVECLIM 1.2 | all components active | ~4 | 100 yr in ~4 h CPU on one Xeon, Goosse et al. (2010) |
| this stack | ExoPlaSim T42 L10, 16 ranks, FP64 | ~76 | 77-93 s per 180.7-day orbit, measured here |
| BIG-MITgcm | SPEEDY L5 + MITgcm L25, CS32 | ~300 | stated directly, Moinat et al. (2026) |

UVic ESCM 2.10 states 4.6 to 11.5 hours per 100 years on a desktop without
naming a core count, so it cannot be placed in this column; on wall clock it
sits between LOVECLIM and this stack. ExoCAM's ~35,000 core-hours is per case
rather than per unit model time and does not convert.

Two readings follow, and they point in opposite directions.

**The EMIC advantage is structural and already spent.** CLIMBER-X is about
twenty times cheaper per simulated year than this stack, and the reason is that
it does not solve the primitive equations: a 2.5-D statistical-dynamical
atmosphere on a 5 degree grid parameterises the circulation a spectral core
computes. For Earth that is a good trade. Here the circulation is one of the
things being asked about, under 30 hour rotation and 32 degree obliquity, so it
cannot be a parameterisation fitted to Earth's.

**A three-dimensional ocean costs about four times.** That is the most useful
single number in the table for OCN-3. BIG-MITgcm carries 25 levels of
primitive-equation MITgcm at an atmospheric resolution comparable to T42 for
about four times this stack's cost per simulated year. cGENIE at 36 x 36 x 16
frictional-geostrophic is a far lighter object than that, so four times is a
pessimistic bound on the ocean plan rather than an estimate of it.

## 3. PALADYN and its neighbours in CLIMBER-X

PALADYN (Willeit and Ganopolski, 2016, `10.5194/gmd-9-3817-2016`) is CLIMBER-X's
land surface model. `references/INDEX.md` records it as read. It matters here
because `references/INDEX.md` already holds JULES as the coupled land-surface
precedent for LSHY, and JULES is CMIP-class; PALADYN is the same object built
for an EMIC, so it is the precedent at a throughput this project shares.

Its shape: nine surface types, five soil layers to 3.9 m with a 20 cm top layer,
one snow layer with prognostic temperature and density, 1-D heat diffusion with
explicit snowmelt and soil phase change, and vertical Darcy flow on
Clapp-Hornberger properties with a free-drainage lower boundary. Physics and
photosynthesis integrate implicitly at a 1 day step, vegetation and soil carbon
at 1 month.

### 3a. The saturated fraction is the other answer to GW-6

PALADYN does not resolve the valley-to-ridge water table gradient either. It
obtains the saturated area statistically, under TOPMODEL (Beven and Kirkby,
1979) as implemented by Niu et al. (2005):

    f_sat = f_sat_max * exp(-f_grad * z_grad)

`f_sat_max` comes from a compound topographic index; PALADYN takes it from
ETOPO1 through Stocker et al. (2014). Wetland fraction is set equal to `f_sat`
wherever the surface is snow free. Infiltration capacity is
`k_sat * (1 - f_sat)`, and runoff is saturation excess over `f_sat` plus
infiltration excess over the rest.

Two things about this are worth stating separately.

**It changes the quantity being asked for.** GW-6 declares the valley-to-ridge
gradient as sub-grid and measured it as unrecoverable at 15.19 km: the model's
water table is a function of recharge at Spearman -0.977 where the real one is
not at -0.156, and GW-23 sized the discarded relief at about three times the
signal being predicted. TOPMODEL carries that terrain control statistically
instead of resolving it, and it delivers a per-cell saturated FRACTION rather
than a per-cell DEPTH. A fraction is what WET-2's wetness classification,
SURF-7's discharge mask and LSHY-6's hydrologic mosaic actually consume. The
depth is the quantity that failed its bar; the fraction has not been tested here.

**The index is computable from the mesh this project already has.** PALADYN
needs ETOPO1 because CLIMBER-X has no sub-grid terrain. Here the 15.19 km region
mesh under a roughly 300 km T42 cell is exactly the population a compound
topographic index is computed over, and GRID-2 is already chartered to persist
the sub-grid distribution. A compound topographic index is not the same
quantity as GRID-2's hypsometry, since it needs upslope contributing area and
local slope rather than an elevation distribution, but both come off the same
mesh and the drainage side already exists in `hydrography`.

Two caveats, neither fatal. `f_sat` still keys on a grid-cell mean water table,
so it inherits whatever bias that mean carries. And PALADYN estimates its mean
water table directly from column volumetric water content rather than from a
lateral solve, which is a much cheaper object than this project's
Dupuit-Forchheimer solver and sidesteps the thing GW-3 and GW-22 found hard;
whether the two should be reconciled or kept as separate estimators is open.

### 3b. Snow compaction carries gravity, and this model has no compaction

PALADYN's snow density is prognostic. Fresh snow density is temperature
dependent after Anderson (1976), and self-loading compaction follows Kojima
(1967) as implemented by Pitman et al. (1991):

    d(rho_sn)/dt = 0.5 * g * rho_sn * w_sn / eta + P_s * (rho_fresh - rho_sn) / w_sn
    eta = eta_0 * exp(k_T * (T_0 - T_sn) + k_rho * rho_sn)

Metamorphism and the effect of melting are neglected.

The compaction term is linear in `g`. At 12.81 against 9.81 m/s2 it runs 1.306
times Earth's for the same load and viscosity, so this planet's simulated
snowpack densifies faster and sits thinner for the same water equivalent.

ExoPlaSim has no compaction at all. Snow density is the constant 330 kg/m3,
declared twice: `rhosnow` at `landmod.f90:73` and `CRHOSN` as a hardcoded
`parameter` at `icemod.f90:24`. Both are used the same way, as
`dsnowz * 1000 / rhosnow`, to convert water equivalent to physical thickness.
That thickness sets the snow thermal layer thicknesses `zsntop` and `zsnowz` in
`landmod`, and in `icemod` it sets the combined ice-plus-snow threshold at
`icemod.f90:591`. So the constant reaches surface albedo through snow cover
depth, the ground heat flux through layer thickness, and the sea ice diagnostic.

330 kg/m3 is a reasonable Earth seasonal-snowpack mean. It is an inherited Earth
constant of exactly the class `notes/audits/inherited-earth-constants.md`
catalogues, it is duplicated with one copy hardcoded beyond namelist reach, and
the mechanism that would make it differ here is one this model does not contain.
The direction is known and one-signed; the magnitude is not, because `eta`
depends on temperature and load and this world's snow temperatures are not
Earth's.

### 3c. A third acceleration, with its invalidity declared

PALADYN can run its vegetation and carbon modules at an artificial internal time
step of 1000 years or more. This is possible only because those components are
fully implicit. In that mode they are called once at the end of each simulated
year, on annually cumulated NPP and litterfall and annual mean decomposition
rates, and the model reaches near equilibrium in about 100 simulated years
instead of the 10,000 the natural timescales would require. 100 years is also
what the physical state needs, permafrost above all.

The paper then states what the mode cannot do: it does not apply to processes
intrinsically out of equilibrium, naming peatlands and inert permafrost carbon,
and getting those to a present state requires a transient run over at least one
glacial cycle.

That declaration is the part worth copying. This tree now has three devices in
the same class and one of them is unrun: `NCVEG` for vegetation carbon,
PALADYN's equilibrium spinup mode, and ExoPlaSim's own `newsnow` for land ice,
whose feedstock every glacier-enabled run here has written and nothing has read
(`notes/audits/dormant-exoplasim-modules.md` finding 5). CLIM-53 asks for a
verdict on the third. PALADYN is the worked example of the verdict's shape: not
whether the accelerator is correct, but what it is valid for and what has to be
integrated honestly instead.

### 3d. The code is reachable, and the licence is clean

CLIMBER-X is public at `cxesmc/climber-x` under GPL-3.0, with PALADYN in
`climber-x/src/lnd/`. The tree has moved past the 2016 paper: lakes are implemented where
the paper calls them a placeholder, and there is now `weathering.f90`,
`dust_emis.f90`, `n2o_emis.f90` and water isotopes.

ExoPlaSim is GPL version 2 "or, at your option, any later version"
(`vendor/exoplasim/LICENSE.TXT` line 298). A GPL-3 combination is therefore
permitted, with the combined work under GPL-3. That answers the licence question
before it is asked, and it matters because LSHY-3's land column would live
inside `landmod.f90`.

The modules that map onto open rows are `soil_hydro.f90` and `surface_hydro.f90`
for LSHY-3, `soil_temp.f90` for LSHY-5, and `water_check.f90`, which is a water
conservation checker of roughly the shape LSHY-7 describes.

### 3e. What is not worth taking, and what is a floor rather than an upgrade

Two of these are refusals and three are something else, and the distinction
matters more than the list.

**Not worth taking: the vegetation.** TRIFFID with five plant functional types,
against LPJ-GUESS CNP. Not close, and in the wrong direction.

**Not worth taking: the uniform hydraulic defaults.** `k_sat`, `psi_sat` and the
Clapp-Hornberger `b` are global uniform values by default, with a texture and
organic-matter formulation available in an appendix. This project has a
lithology-derived soil, so only the second path is relevant, and LSHY-1 already
owns the property contract.

**A floor, not an upgrade: peat structure, peatland area and the methane
fraction.** WET-4, WET-6 and WET-9 charter this ground to a standard PALADYN
does not meet: vertically resolved substrate, temperature, saturation and redox
state, production and oxidation exposed independently, and stock-conserving
transition semantics. Measured against that specification PALADYN is a step
down, and adopting it in place of the specification would be a step down from a
standard already set.

That comparison is against a TARGET, though, and the target is unproven. WET-1
through WET-11 specify what an accepted wetland, peat and methane component
would be; nothing has yet shown the chain is implementable within this project,
and it rests on SDEC-3's vertical coordinate, PLHY-4 and PLHY-5, PCAR-5's trait
registry and BVOC-6's oxidant bracket, all of which are wholly open. A standard
that does not ship is worth less than a design that does. PALADYN is a design
that works, is published, and runs at this project's throughput, so it is worth
recording as the floor each part can fall back to rather than as a rejected
option.

The three parts are not equally viable as a floor, and this is the substance:

- **Extent is genuinely supplied.** Wetland fraction is the TOPMODEL saturated
  fraction wherever the surface is snow free, and potential peatland is the
  fraction wet for at least three months of the year. That is a geographic
  quantity computed from terrain and column water, and section 3a's index is its
  input. It forfeits nothing WET-2 needs except the mutual exclusivity WET-2
  asks for, which is a bookkeeping requirement rather than a physical one.
- **Peat stock is supplied, at a cost that is nameable.** Acrotelm and catotelm
  confined to the top soil layer, transfer at a critical acrotelm carbon of
  5 kgC/m2 after Wania et al. (2009), catotelm shifted to lower layers as
  density fills one, peat carbon vertical diffusivity set to zero, and areal
  change after Stocker et al. (2014) limited to 1 percent a year with a minimum
  fraction seeding every cell. What it forfeits is WET-4's vertical resolution,
  and with it any depth-resolved redox claim.
- **Methane is a bracket and never a measurement.** A constant fraction of
  heterotrophic respiration wherever respiration is anaerobic, one fraction per
  surface type. It cannot separate production from oxidation, which is what
  WET-6 exists to expose, and it contains no transport at all, so WET-7's
  diffusion, plant and ebullition pathways collapse into the constant. It yields
  a surface FLUX under a transplanted-Earth assumption. It does not yield an
  abundance, because an abundance needs the oxidant and lifetime calculation
  WET-10 owns and this form does not contain, so it can never move
  `atmosphere.pCH4_bar`.

WET-12 carries the declaration and the condition for taking each part.

### 3f. Three more CLIMBER-X modules, and one of them fills the hole section 4 leaves

PALADYN is one directory of `cxesmc/climber-x`. Three of its neighbours matter
here, and the licence position of section 3d covers all of them. Paths below
are in THAT repository and not in this one, which is why they are written with
the `climber-x/` prefix: a bare `src/` in this tree reads as `docs/src/`, and
`source/` is the Orogen exports directory that rule 7 governs.

**`climber-x/src/smb` is a surface energy and mass balance model, and it is the answer to
what section 4 refuses.** MITgcmIS's Positive Degree Day scheme is rejected
below because it drives ablation from 2 m air temperature alone. SEMI is the
opposite kind of object. Its interface takes surface albedo and downward
shortwave in FOUR components each, visible and near-infrared by direct and
diffuse, together with the derivatives of downward shortwave with respect to
albedo, plus cloud cover, downward longwave, a lapse rate, a cosine of zenith
angle and dust.

That band structure is this project's own. The visible and near-infrared split
is the same decomposition `lib/stellar.py` computes at 0.75 um and that
`world_state.json` records as `flux_fraction_band1`, so the per-star reweighting
already done maps onto SEMI's inputs instead of having to be invented. A scheme
built for a two-band surface is the one kind of mass balance model a non-solar
host does not immediately break.

It also takes sub-grid surface elevation, an elevation standard deviation,
surface slopes and elevation classes, and calls `downscaling_mod` for radiation,
precipitation and wind. So it runs the surface balance on a DOWNSCALED field
rather than on the cell mean. That is a third independent answer to the problem
PHYS-13 and GRID-2 describe, and unlike the freezing-height criterion, which is
a diagnostic, this one is a scheme that consumes the sub-grid distribution as an
input. `smb_surface_par.f90` and `snow.f90` carry the surface side, and
`snow_par` exposes `lsnow_dust`, `w_snow_dust` and `dust_con_scale`, so dust
darkening of snow is a namelist switch rather than a gap. That is DUST-14,
already solved in this model class.

Two things in it do not transfer. Every constant is a PHYS-class inheritance,
and `smb_bias_corr.f90` is a bias correction against Earth observations, which
has no meaning on a world with none.

**`climber-x/src/ch4` is a reduced atmospheric methane model, and WET-10 has a floor after
all.** `ch4_model.f90` is 12 kB and is a partitioned-lifetime box: separate
tropospheric OH, chlorine, soil and stratospheric sink timescales, an OH
temperature sensitivity, and OH sensitivities to CO, NOx and VOC precursors,
with emissions converted to abundance in ppb.

The structure is what transfers and the calibration is not. The precursor
sensitivities are fitted to anthropogenic Earth emissions read from a file, the
chlorine term is Earth stratospheric chemistry, and the OH field itself is a
photochemical product of the host's ultraviolet. For this star that last one is
not an open question: `config/planet.yaml` records that Rugheimer et al. (2013)
computed the photochemical steady state at this effective temperature, which is
the same calculation an OH field would come out of. So WET-10's oxidant and
lifetime half has a reduced form available in shape, with the coefficients as
the work.

**`climber-x/src/sic` is dynamic-thermodynamic sea ice**, `sic_dyn.f90` and
`transport_sic.f90`, against this project's thermodynamic-only scheme.
BIG-MITgcm names excessive ice as the consequence of the configuration this
stack also runs, so the term is not free. It is nevertheless BLOCKED rather than
available: sea ice dynamics needs wind stress and ocean currents, and this ocean
is a slab with horizontal transport deliberately off. It becomes a question when
the ocean plan lands and not before.


## 4. BIG-MITgcm is the nearest peer, and two of its details are warnings

BIG-MITgcm (Moinat et al., 2026, `10.5194/gmd-19-4357-2026`,
`references/gmd-19-4357-2026.pdf`) is the only published model that made the
same architectural bet as this project. SPEEDY physics on the MITgcm core at
CS32, about 2.8 degrees, with a 25-level ocean, Winton thermodynamic sea ice, a
two-layer land model, BIOME4, pysheds runoff routing and a new shallow-ice
module. Coupling is asynchronous and offline: run the fast system to a steady
state, update the slow components, repeat until convergence. That is loop A, B
and C with different components in the slots.

Its convergence predicate is declared: fewer than 10 percent of land points
changing biome or total ice volume, and a surface energy imbalance below
0.2 W/m2. Five iterations were needed, each requiring several thousand simulated
years to reach that steady state.

**Its ice sheet module is nearly free, and it splits in two.** MITgcmIS is
shallow-ice with Glen's law at n = 3, neglecting basal sliding, calving and
basal melt as unresolvable at that resolution, and it costs about 1 CPU-hour per
40,000 years. It also carries LLRA isostatic adjustment and lapse-rate,
freshwater and sea-level corrections. The code is on Zenodo at
`10.5281/zenodo.18723952`. So ice sheet flow is a small component rather than a
large one, which is the useful correction to make.

The FLOW half is worth taking and would not be ported. Equation 4 is

    dH/dt = div(D grad H) + div(D grad z_B) + A_dot
    D = 2a/(n+2) * (rho_i g)^n * |grad z_S|^(n-1) * H^(n+2)

which is a few lines of algebra, and `(rho_i g)^n` is exactly the term GRAV-6
discusses at 2.23x here, arriving correctly rather than bolted on.

There are two routes to it and this note does not choose between them. MITgcmIS
is Python on MITgcm's cubed sphere while this project's grids are Gaussian, and
CLAUDE.md rule 3 exists because that class of crossing has silently matched zero
cells three times, which argues for implementing on the native mesh from the
paper instead. Against that, a doubly nonlinear diffusion has direct local
precedent for failing to converge on this mesh: GW-9's exponential
transmissivity failed under both Picard and Kirchhoff and was abandoned. Which
route is cheaper has not been measured and is CLIM-62's; nothing here should be
read as having settled it.

The MASS BALANCE half is where the difficulty lives, and theirs is weaker than
what this project already plans. Its only two inputs are 2 m air temperature for
ablation and snow precipitation for accumulation, taken as daily output over 30
years and averaged per day and per cell from a steady state.

- **Ablation is temperature-only.** The Positive Degree Day method after Tsai
  and Ruan (2018) is better than a bare degree-day factor, since a percolation
  layer of thickness `H_p` gives a semi-physical melt delay, but no shortwave,
  albedo or spectrum enters it. Under a K2.5V host that is the term this world
  differs on most over ice: `config/planet.yaml` records that ExoPlaSim's own
  `k2.dat` made snow and ice 0.10 to 0.17 too dark in every run before the
  measured spectrum replaced it. A temperature-only ablation scheme is
  structurally blind to the quantity that correction exists to get right.
- **Accumulation assumes away GRAV-8.** The paper states that the MITgcm land
  module has no process that densifies snow into glacial ice, so densification
  is assumed instantaneous and snow precipitation is divided directly by
  `rho_i = 920 kg/m3`. That is the same missing physics GRAV-8 opens, handled by
  assuming it does not matter.

**And the ordering binds regardless of either.** `glac` is identically zero
across all twelve bins and all 8,192 cells, and PHYS-13 establishes that this is
a grid artifact rather than a fact about the planet: corrected to each mesh
region's elevation, land below freezing in the warmest month goes from 0.002 to
1.657 percent, because the model evaluates its own high ground about 7.8 K too
warm on average and 21.9 K in the top tenth. A temperature-driven ablation
scheme is the most sensitive possible consumer of exactly that bias, so run
today it would return zero ice, confidently, for the wrong reason. The sequence
that makes the question answerable is PHYS-13's ice mask, then CLIM-53's
accelerator verdict, and only then whether flow is needed at all: `newsnow`
grows ice without moving it, which is the cheaper test of whether ice persists
for a year. The one standing claim that will eventually need flow is
`config/planet.yaml`'s, that the long stellar cycle component is slow enough
that the simulation's glaciers equilibrate rather than merely breathing.

**Its cloud albedo is tuned on latitude for agreement.** The description is
explicit: cloud albedo depends on latitude in order to reduce net solar
radiation at high latitudes and therefore to agree better with observational
data. That is a fitted knob on a radiative quantity in a peer-reviewed model
description, and it is what `docs/src/practice/failure-modes.md` class 16
forbids by name. Anything borrowed from a SPEEDY-derived scheme carries it.

**Its bathymetry preparation deletes closed water.** Isolated oceanic points and
lakes are removed irrespective of size to avoid numerical instability, and every
point shallower than 20 m is set to 20 m. On a world where about three quarters
of land drainage is endorheic, a preparation step that silently removes enclosed
water is not a detail. OCN-11 already requires bathymetric smoothing and
connectivity changes to be declared explicitly; this is that requirement
justified by an observed case.

Two of its stated limitations validate diagnoses already made here. Its coarse
grid underestimates Greenland and Antarctic elevation and runs those regions too
warm, which is GRID-2 and PHYS-13's sub-grid orography problem found
independently at a comparable resolution. And its thermodynamic-only sea ice
produces excessive ice, which is the configuration this stack also runs.

Its atmosphere resolves longwave in four bands against this model's Sasamori
broadband plus the one added trace-gas band. On radiation this peer is ahead,
and it points the same way `exoplasim/notes/corrk-cross-check.md` does.

## 5. CliMA is not a competitor and has one thing worth reading

The Climate Modeling Alliance is writing a full Earth system model in Julia for
GPUs: `ClimaAtmos.jl`, `ClimaOcean.jl` on Oceananigans, `ClimaLand.jl`,
`ClimaSeaIce` and `ClimaCoupler.jl`, with correlated-k radiation in
`RRTMGP.jl`. Every one of those repositories had commits within the month of
this survey. It is CMIP-class, so its throughput is three to four orders below
the band this project works in, and nothing in that regime can carry a terrain
loop.

What is worth reading is `ClimaParams.jl`. It carries a block headed *Planetary
and Orbital Parameters (Earth Defaults)*: `planet_radius`,
`gravitational_acceleration`, `angular_velocity_planet_rotation`,
`length_orbit_semi_major`, `total_solar_irradiance`, `orbit_obliquity_at_epoch`,
`orbit_eccentricity_at_epoch` and a longitude of perihelion. That is
structurally `config/planet.yaml`, including the decision to name Earth's values
as defaults rather than compile them in. No EMIC has anything like it, and it is
an externally maintained enumeration of exactly the quantities OCN-12's
implicit-Earth audit is looking for.

`RRTMGP.jl` also treats the stellar source function as data rather than a
constant, though its correlated-k tables are trained on Earth-like gas mixtures
over Earth-like abundance ranges.

## 6. Candidates this survey closes

**A resolved-tier ocean in the model tree.** `lsgmod.f90` ships with ExoPlaSim
and is doubly off: `compile.sh:305` exports `OCEANCOUP=cpl_stub`, so every
coupling call reaches an empty return, and the full LSG source in `lsg/src/` is
never built. OCN-3 removed the resolved tier on spin-up cost and a coarse
geostrophic ocean does not escape that argument.
`notes/audits/dormant-exoplasim-modules.md` records this so it is not
rediscovered as an available option; nothing in this survey reopens it.

**A GPU port of the spectral transform.** Priced and parked under CLIM-58 on
this machine's arithmetic. CliMA is the existence proof of the other answer,
which is to design for GPUs from the first line and treat mixed precision as a
modelling decision. That is a rewrite rather than an optimisation and is not on
this project's table.

**Vendoring PALADYN whole.** A fourth subtree is a maintenance commitment, and
the Earth configuration inside it is the OCN-12 problem in a new place: uniform
hydraulic defaults, an ETOPO1 topographic index, calibrated methane fractions
per surface type. The formulations in section 3 are the valuable part and they
are all in the paper.
