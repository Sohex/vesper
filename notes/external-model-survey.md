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

## 2. Cost, and why the obvious comparison does not hold

**The table below compares things that are not the same object, and every
conclusion originally drawn from it was wrong in the same direction.** It is
kept, with the error stated, because the error is instructive and because the
numbers themselves are reported correctly by their sources.

| model | what the figure covers | CPU-h per 100 sim-yr | source |
| --- | --- | ---: | --- |
| CLIMBER-X | atmosphere, 3-D ocean, sea ice AND the PALADYN land surface | ~3.8 | 10,000 sim-yr/day on 16 CPUs, Willeit et al. (2022) |
| LOVECLIM 1.2 | ECBilt, the CLIO ocean and VECODE vegetation, all active | ~4 | 100 yr in ~4 h CPU on one Xeon, Goosse et al. (2010) |
| this stack | ExoPlaSim ALONE, on a slab ocean | ~76 | 77-93 s per 180.7-day orbit, measured here |
| BIG-MITgcm | SPEEDY L5, a 25-level ocean, sea ice and a 2-layer land model | ~300 | stated directly, Moinat et al. (2026) |

UVic ESCM 2.10 states 4.6 to 11.5 hours per 100 years on a desktop without a
core count and cannot be placed in this column. ExoCAM's ~35,000 core-hours is
per case rather than per unit model time and does not convert.

### 2a. The error

Row three is one component. Rows one, two and four are coupled systems. Reading
down the column and concluding that this stack "lands inside the EMIC band"
compares an atmosphere against suites that carry an ocean, sea ice, a land
surface and in two cases vegetation.

The comparison is too generous in both directions at once. It flatters this
stack, because a Vesper commissioning is not one ExoPlaSim run: it is a
bootstrap and a baseline, hydrography, pedology, an hours-class MPI LPJ-GUESS
run, the aeolian components and a carve verdict, and an ITERATION adds an Orogen
generation in front of all of it. `docs/src/pipeline/costs.md` already prices the
work in those units and this section did not. And it flatters CLIMBER-X and
LOVECLIM less than it looks, since their component counts are what their small
numbers are buying.

**The 4x reading was the worst of it.** "A three-dimensional ocean costs about
four times" came from 300 against 76, and that ratio is a whole coupled system
against a bare atmosphere, not the price of an ocean. It is also confounded on
the atmosphere side: SPEEDY at five levels on CS32 is a much cheaper atmosphere
than ExoPlaSim at T42 with ten. Two unlike things differ, and the difference was
attributed to the one term that happened to be interesting.

What partially survives is narrower and worth keeping. BIG-MITgcm's own
breakdown puts BIOME4 and pysheds under five minutes and MITgcmIS at about
1 CPU-hour per 40,000 years, so its 300 is almost entirely the coupled fast
system rather than the slow components. That makes it an atmosphere-plus-ocean
figure against an atmosphere-plus-slab figure, which is a real comparison with
one confound rather than none.

### 2b. The unit does not transfer, which is the deeper problem

CPU-hours per 100 simulated years assumes every component advances on the same
clock. This pipeline is asynchronous by construction: LPJ-GUESS runs once per
pass of loop B, the carve verdict once per pass of loop A, the aerosols once per
climatology. There is no rate at which the biosphere runs "per simulated year",
so the unit is not merely inconvenient here, it is undefined.

The units this project already uses are the right ones and they have no
counterpart in the table: cost per commissioning, and cost per iteration.
CLIMBER-X and LOVECLIM have no equivalent because nothing in them regenerates
the terrain.

So the honest position is that this stack's throughput has NOT been compared
with the EMIC family, and cannot be by that column. What has been measured is
that one ExoPlaSim orbit costs 77 to 93 seconds at T42 on 16 ranks, which is a
fact about the atmosphere and is used correctly elsewhere in this note. OCN-3
still needs cGENIE's cost, and SPAT-11 still needs the per-rung cost; neither is
answered by the row above.

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

CLIMBER-X is public at `cxesmc/climber-x` under GPL-3.0-or-later, and a pinned
extraction is held locally at `references/climber-x/`, with PALADYN in
`references/climber-x/src/lnd/`. `references/INDEX.md` records the revision. The tree has moved past the 2016 paper: lakes are implemented where
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
here, and the licence position of section 3d covers all of them. Paths below are on
disk under `references/climber-x/` and are written in full: a bare `src/` in
this tree reads as `docs/src/`, and `source/` is the Orogen exports directory
that rule 7 governs.

**`references/climber-x/src/smb` is a surface energy and mass balance model, and it is the answer to
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

**`references/climber-x/src/ch4` is a reduced atmospheric methane model, and WET-10 has a floor after
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

**`references/climber-x/src/sic` is dynamic-thermodynamic sea ice**, `sic_dyn.f90` and
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
freshwater and sea-level corrections. The code is on Zenodo at `10.5281/zenodo.18723952` and is held locally at
`references/big-mitgcm/`, where `MITgcmIS.py` is the whole ice sheet model in
557 lines; its gigabyte of simulation output was not fetched. So ice sheet flow is a small component rather than a
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
forbids by name.

It is THEIR modification and not SPEEDY's, which is worth stating because the
distinction decides whether it travels. Stock SPEEDY applies two constants,
`albcl = 0.43` and `albcls = 0.50` in `shortwave_radiation.f90`, uniformly and
with no latitude term; the latitude dependence is BIG-MITgcm's, after Ragon et
al. (2022). So the knob is local to their configuration, and the four longwave
bands section 8a wants are not carrying it.

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

**Surveying the other EMICs' source.** LOVECLIM, UVic ESCM and Bern3D were
compared on cost and component list and their source was not opened. That is
deliberate rather than unfinished: the value in CLIMBER-X is concentrated in its
LAND SURFACE, which is the one place an EMIC could be ahead of this project, and
those three are weaker there than CLIMBER-X by construction. LOVECLIM's
vegetation is VECODE at two plant functional types, and UVic and Bern3D pair
energy-moisture-balance atmospheres with land schemes simpler still. Where they
lead is ocean biogeochemistry, and OCN-3 and OCN-4 already own that question
with cGENIE and ECOGEM named. Reopening this needs a reason from a specific
open row, not a second sweep.

**Vendoring PALADYN whole.** A fourth subtree is a maintenance commitment, and
the Earth configuration inside it is the OCN-12 problem in a new place: uniform
hydraulic defaults, an ETOPO1 topographic index, calibrated methane fractions
per surface type. The formulations in section 3 are the valuable part and they
are all in the paper.


## 7. REF-10: the rest of the CLIMBER-X tree, triaged

*Read 2026-08-22 against `references/climber-x/` at the pinned revision. The
rule was declared in REF-10 before the reading: a module is a HIT when it
implements a mechanism an open row already names, and is recorded and dropped
otherwise. Depth varies deliberately. Where a directory listing was enough to
apply the rule, that is all it got, and this says which.*

### 7a. Hits, keyed to the rows that were waiting

**`src/ocn` -> OCN-19, OCN-20, OCN-5, OCN-12.** The ocean shares GOLDSTEIN's
lineage with cGENIE, and it is not a family resemblance: `invert.f90:70` orders
the streamfunction points as `k = i + j*n`, which is what OCN-20 quotes from
cGENIE's `invert.f:32`, and `ubarsolv`, `island`, `matinv` and `jbar` are all
still called from `momentum.f90`.

What has changed matters to those two rows. `ubarsolv` has been reworked for
CONTIGUOUS BAND STORAGE, with the LU factors repacked as `Lband`, `Uband` and
`Udiag` built once in `momentum`, and a unit-stride dot product over the upper
band in back substitution. It is still two triangular sweeps, so the `r^3` per
timestep OCN-20 targets is untouched. **The memory-layout half of that
optimisation exists in a GPL-3 sibling and the complexity half does not**, which
splits OCN-19 and OCN-20 more cleanly than they were split when written.

`free_surface.f90` is called from `ocn_model.f90:533` but diagnoses sea surface
height from density; it does not replace the rigid-lid machinery, and reading
the filename as though it did would have been wrong. `restore_salinity.f90`,
`flux_adj.f90` and `hosing.f90` are the named controls OCN-5 requires be kept
explicit rather than physical, and `eos.f90` is an OCN-12 item.

**`src/geo/hypso_topo.f90` -> GRID-2.** A worked implementation of the shape
GRID-2 declares: a high-resolution bed elevation binned into per-coarse-cell
AREA FRACTIONS, 10 m bins to 6500 m, with an update cadence and a finer
refinement band over a depth range of interest. It is oriented at ocean depth
for sediment and coral rather than at land relief, so the binning is not
transferable and the shape is.

**`src/geo` hydrology -> HYD, OCN-11.** `lakes.f90` is a dynamic lake model,
`topo_fill.f90` fills topography, and `drainage_basins.f90`, `runoff_routing.f90`
and `fix_runoff.f90` sit beside them. Depression filling is the algorithmic
counterpart to the carve list. `coast_cells.f90`, `connect_ocn.f90` and
`fill_ocean.f90` are ocean connectivity and are OCN-11's contract implemented,
which is worth reading beside BIG-MITgcm's delete-the-lakes step rather than
instead of it.

**`src/main/constants.f90` -> OCN-12, and it is a negative hit of the most
useful kind.** `R_earth`, `omega`, `fcoriolis = 2*omega`, `g = 9.81` and the
WGS84 ellipsoid are Fortran `parameter` constants, and so is
`frac_vu = 0.45`, commented "fraction of solar spectrum in visible and
ultraviolet". That is the Sun's band split hardcoded at compile time, and it is
the exact opposite of `ClimaParams.jl`'s declared planetary block in section 5.

It also PRICES the standing advice to take formulations and not the tree. Every
module recommended in section 3 imports Earth constants from this file except
one: `surface_hydro.f90`, `soil_temp.f90` and `surface_par_lnd.f90` all do, as
do every ocean module, all of `geo`, and `sico_params.f90`. `soil_hydro.f90`
does not. So severing that dependency is part of the cost of any take, and this
file is its audit list.

**`src/ice_sico` -> CLIM-62.** A third route, and the heaviest: SICOPOLIS with
full thermomechanics, a temperature-dependent rate factor, an enhancement
factor and finite-viscosity regularisation on Glen's law at n = 3. That is
exactly what BIG-MITgcm omits and what GRAV-6 notes is sometimes included.
Recorded as an option with its weight named; nothing here recommends it.

**`src/bnd/fake_*.f90` -> OCN-10, EFOR-1.** The offline-driver pattern, worked:
`fake_lnd.f90` reads runoff and discharge from one netCDF file, precomputes
monthly-to-daily interpolation weights, and presents the same derived type the
real component would, so a component can be swapped for a file without its
consumers knowing. There are eight of these, one per component.

**`src/main/coupler.f90` -> a principle rather than a row.** Every exchange goes
through one COMMON GRID, `cmn_to_atm` and `atm_to_cmn` and their siblings, never
component to component pairwise. That is CLAUDE.md rule 3's positive form as
architecture. It also confirms section 3f's reading of SEMI from the other side:
`alb_vis_dir_ice_semi`, `alb_vis_dif_ice_semi`, `alb_nir_dir_ice_semi` and
`alb_nir_dif_ice_semi` are carried on the common grid, so the four-component
albedo is the coupled interface and not an internal convenience.

### 7b. Read and dropped

**`src/atm/lwr.f90` and `swr.f90` -> CLIM-61, and the answer is no.** This was
expected to be a middle rung between the broadband scheme and correlated-k. It
is not. The shortwave is TWO bands, visible-plus-ultraviolet against infrared,
split at the hardcoded `frac_vu`, which is the same structure PlaSim already
has. The longwave is not band-resolved at all: it is a fitted parameterisation
over gas concentrations producing a CO2 equivalent, with coefficients like
`ak_o3 = 0.6` carrying the comment "from tuning of total LW contribution by
O3". So SPEEDY's four longwave bands remain the only cheaper rung identified,
and CLIM-61's candidate list does not grow.

**The rest of `src/atm`** is the statistical-dynamical core, `adifa`, `crisa`,
`u2d`, `u3d`, `wvel`, `slp`, `synop`, `vesta`, `diffuse_impl`. This project runs
a spectral primitive-equation core and is not replacing it. `feedbacks.f90` and
`rad_kernels.f90` are radiative-kernel feedback decomposition with no open row
waiting; CLIM-1 closed the energy decomposition question at this resolution.

**`src/co2` and `src/n2o`** are the same box-model shape as `src/ch4`. Neither
gets a row: CO2 is prescribed by decision under OCN-16 and N2O is prescribed
with its abundance measured from Rugheimer, so unlike methane there is no row
waiting for a reduced form to fill.

**`src/bnd`, the rest.** `luc.f90` is land use change, `cfc.f90`, `d13c_atm.f90`
and `D14c_atm.f90` are Earth isotope and halocarbon forcings. `insolation.f90`,
`solar.f90` and `o3.f90` are owned here already by `lib/orbit.py`,
`lib/stellar.py` and a determined `ozone_scale`, with SPEC closed at 0 of 5.

**`src/bmb`** is basal mass balance and is downstream of ice existing at all.
**`src/utils`** is a tridiagonal solver, a filter, a precision module and a
hysteresis helper. **`src/lndvc`**, **`src/ice`** and **`src/bgc-dummy`** are
coupling shims and disabled stubs; `bgc-dummy` is four dotfiles.

### 7c. What this triage did not do

It read for mechanisms against open rows and nothing else. It did not evaluate
correctness, did not build anything, and did not compare numerical results,
because the tree is comparison material rather than a dependency. Six
directories got a listing and a judgement rather than a read, named in 7b, on
the grounds that no open row named a mechanism they contain: a listing is enough
to apply the declared rule and not enough to claim anything else about them.


## 8. SPEEDY and SOCRATES, read for CLIM-61

*Read 2026-08-22 from `references/speedy/` and `references/socrates/` at the
pinned revisions. CLIM-61's costing is ordered after the SHTns forward direction
lands, because both touch the same worktree. Reading is not, and what follows
changes what the row is about.*

### 8a. SPEEDY: the cheap rung is real, and it is only the longwave half

Its SHORTWAVE is no better than what this model already has.
`shortwave_radiation.f90` splits incoming flux at `fband2 = 0.05` and applies
fitted absorptivities per constituent, `absdry`, `absaer`, `abswv1`, `abswv2`,
`abscl1`, `abscl2` and an ozone fraction `epssw`. That is the same class of
object as Lacis and Hansen: two bands and fitted absorptances. Adopting it would
buy nothing and would need the same per-star reweighting `config/planet.yaml`
already carries, because `fband2` is a solar split.

Its LONGWAVE is four bands and that is the part worth having.
`longwave_radiation.f90` declares `nband = 4` and `mod_radcon.f90` carries
`fband(100:400,4)`, the energy fraction emitted in each band as a function of
temperature, against this model's single Sasamori broadband plus the one band
added for CH4 and N2O.

**The longwave partition is star-independent, and that is the whole reason it
transfers.** It divides the PLANET's own emission, so it is a property of the
Planck function at terrestrial temperatures rather than of the host. Nothing in
it needs the treatment every shortwave term here has had.

Two caveats, both fixable and both worth knowing before the row is costed. The
partition is not a Planck integral but a fit:

    fband(T,2) = (0.148 - 3.0e-6*(T-247)^2) * eps1
    fband(T,3) = (0.356 - 5.2e-6*(T-282)^2) * eps1
    fband(T,4) = (0.314 + 1.0e-5*(T-315)^2) * eps1
    fband(T,1) = eps1 - the other three

And it is CLAMPED FLAT outside 200 to 320 K: below 200 K every band takes its
200 K value and above 320 K its 320 K value. Stratospheric temperatures sit
below 200 K routinely, so the clamp is reached in normal operation rather than
at an extreme. Re-deriving those four quadratics from an actual Planck integral
over the band edges is cheap and would remove the clamp; it is the kind of work
`shortwave_band_weights.py` already does for the other half of the spectrum.

### 8b. SOCRATES: the star is a data file, which is the whole argument

The cost first, because it is smaller than the four orders of magnitude quoted
when comparing whole models. In the `ga7` configuration the shortwave file
carries 6 spectral bands, 8 gaseous absorbers and at most 12 k-terms in a band;
the longwave carries 9 bands and 12 absorbers. Against 2 shortwave bands and one
broadband longwave that is roughly an order of magnitude more radiative work,
not four. With `radstep` already the largest single term in the profile at over
12 percent, a tenfold radiation cost is about a twofold whole-model slowdown.
That is a number for CLIM-61 to MEASURE rather than a reason to stop, and the
cost is tunable: `sbin/Ccorr_k` builds spectral files, so dropping absorbers
this world has no use for and reducing k-terms are both available.

**The distribution already supports non-solar hosts and ships the worked
example.** `data/solar/trappist1` is a stellar spectrum, and its header reads
"BT-Settl, teff = 2600 K, logg = 5, meta = 0" -- the same library and the same
three parameters `build_stellar_spectrum.py` uses to produce this star's
spectrum at 4965 K. `examples/trappist1/mk_ga_trappist` is 46 lines. There is
also `examples/mars` and `examples/titan`.

What that script does is the finding. It takes an existing shortwave spectral
file, runs `prep_spec`, points it at the stellar spectrum, and writes a new
spectral file. **The k-distributions are not regenerated.** Gas absorption is a
property of the gas and not of the star; what changes for a different host is
how much incident flux falls in each band, which is a re-weighting of a file
that already exists. And the longwave file is copied into the hybrid unchanged,
which is section 8a's point arriving from the other direction.

**So the argument for CLIM-61 is not only accuracy.** This project carries seven
derived per-star corrections onto a broadband scheme -- `ozone_scale`,
`ozone_uv_weight`, `ozone_visible_weight`, `h2o_sw_weight`, `h2o_sw_level`,
`co2_sw_weight` and `cloud_absorption_scale` -- plus a hand-added longwave band,
and every one of them exists because a fitted broadband absorptance cannot be
re-weighted structurally. You have to re-derive each one. In a band-resolved
scheme the same operation is a tool invocation against a spectrum file, and it
is the operation those seven factors are each a manual instance of.

That reframes the row. The question is not only whether correlated-k is more
accurate at a price. It is whether the seven factors keep needing maintenance
every time the star, the composition or a surface endmember moves, and what
retiring that maintenance is worth against a whole-model slowdown CLIM-61 has
yet to measure.


## 9. The exoplanet models, read without a row to justify it

*Read 2026-08-22 from `references/exort/`, `references/exocam/` and
`references/isca/`. This section exists because the rule REF-10 declared was
wrong, and the correction found something.*

REF-10's rule was that a module is worth reading when an open row already names
its mechanism. That rule cannot find a blind spot, because a blind spot is by
definition a thing no row names. It was also written here rather than being a
standing convention of this project, and then cited as though it were one.

The blind spot in the survey itself was visible once the rule was dropped: three
EARTH models had been read closely, PALADYN, BIG-MITgcm and cGENIE, and no
exoplanet model had been read at all. For a project whose distinguishing feature
is not being Earth, that is the wrong place to have spent all the attention.

### 9a. What the exoplanet models converge on, and it is this project's method

`references/exort/data/solar` ships stellar spectra as BT-Settl models binned
onto ExoRT's band grids, and `tools/makeStellarSpectrum_fromSED.py` bins a raw
SED onto any of them. SOCRATES ships `data/solar/trappist1`, also BT-Settl, with
`teff`, `logg` and metallicity in its header. `build_stellar_spectrum.py` here
takes BT-Settl and this star's effective temperature and metallicity, deriving
surface gravity and refusing a grid point that disagrees.

Three independent codes, the same library and the same three parameters. That
is a convergence rather than a finding, and it is worth recording precisely
because nothing needs to change: the stellar input this project builds is the
input both of those codes would want.

ExoRT also gives CLIM-61 a second cost point with structure SOCRATES does not
have. It ships FOUR band configurations, n28, n42, n68 and n84, chosen by
atmospheric regime. So correlated-k is a ladder rather than a single price, and
`ga7`'s 6 shortwave and 9 longwave bands is one rung of several.

### 9b. Snow and ice albedo were never re-derived, and nothing had noticed

`references/exocam/tools/py_progs/broadband_albedo_calculator.py` integrates a
reflectance spectrum against a stellar spectrum to a broadband albedo. That is
what `analysis/vegetation_albedo.py`, `analysis/playa_albedo.py` and
`analysis/rock_albedo.py` each do by hand, so the tool is a second
implementation to check those three against, which is the kind of test the
conventions ask for and this project cannot otherwise construct.

Looking for the fourth one is what found the gap. There is no
`analysis/snow_albedo.py`. `config/planet.yaml` carries no snow or ice albedo
key. Nothing here writes `dsnowalbmn`, `dsnowalbmx` or `dglacalbmn`, so all
three take ExoPlaSim's namelist defaults while vegetation, playa clastics and
every rock class carry values reweighted to this star.

The argument that those three needed reweighting applies to snow unchanged, and
in the opposite direction. Snow's reflectance is a property of snow, but band 2
runs from 0.75 um into the infrared, snow darkens steeply across it, and this
star puts 61.8 percent of its flux there. A leaf is dark in the visible and
bright in the near infrared, so this star RAISED the vegetation endmember; snow
is the reverse and should fall.

What makes it worth a row rather than a note is leverage. `config/planet.yaml`
already records that the wrong stellar spectrum made snow and ice 0.10 to 0.17
too dark in every run before it was replaced, so the sensitivity is established
and measured. The spectrum was fixed. The endmembers were not. And snow and ice
albedo drive the ice-albedo feedback, which is the dominant feedback at the cold
end, so this reaches the bare-rock endmember arm, PHYS-13's ice mask, the
extreme-cold land fraction the design flux is capped by, and what the stellar
cycle's grand minima are worth. PHYS-14 owns it.

`references/exocam/tools/spectral_albedos` ships `snow100um.txt`,
`Bluemarineice.txt` and a mixture, so the reflectance data is available rather
than to be found.

### 9c. Isca has already wired SOCRATES into a spectral GCM

`src/atmos_param/` is a ladder of interchangeable physics: `two_stream_gray_rad`,
`rrtm_radiation`, `sea_esf_rad`, several convection schemes, and `qflux`, which
is a prescribed ocean heat flux of the kind OCN-2 verifies here. And
`src/atmos_param/socrates`.

**That module is the answer to CLIM-61's interface question**, which is most of
the work in adopting a band-resolved scheme. The whole integration is 3,452
lines, of which `socrates_interface.F90` is 1,688 and `socrates_set_cld.F90` is
503, and the spectral files arrive as namelist strings, `sw_spectral_filename`
and `lw_spectral_filename`, so pointing it at a K dwarf file is configuration
rather than code. Isca is FMS-based and ExoPlaSim is not, so nothing drops in;
what transfers is the SPECIFICATION, because the call signature states exactly
what a host model has to supply:

    temperature, specific humidity, ozone, CO2, surface temperature,
    pressure at full and half levels, height at full and half levels,
    albedo, cosine of zenith angle, Sun-planet distance factor,
    cloud fraction, effective radius, cloud water content

and it returns heating rates, up and down fluxes, clear-sky counterparts, and
optionally a spectral outgoing longwave.

One detail in that list matters for PHYS-14. The albedo argument is a SINGLE
broadband value per column, not the four-component visible and near-infrared by
direct and diffuse that SEMI takes. So a band-resolved atmosphere does not by
itself resolve the surface, and what surface albedo is fed to it remains a
decision rather than something the scheme answers.

### 9d. An exoplanet code with an Earth constant in it, and the check it prompted

`src/shared/constants/constants.F90` puts the planetary constants in a runtime
namelist, `constants_nml`, carrying `radius`, `grav`, `omega`,
`orbital_period`, `rotation_period`, `solar_const`, `pstd`, `rdgas`, `kappa`
and more, with `omega` derived as `2*pi/rotation_period` rather than set
independently. That places Isca with `ClimaParams.jl` and with
`config/planet.yaml`, and leaves CLIMBER-X's compile-time `parameter` block as
the outlier of the four codes read here.

But not entirely. `RADCON`, the conversion from radiative flux divergence to a
heating rate, is a `parameter` built from `EARTH_GRAV`, `EARTH_CP_AIR` and
`SECONDS_PER_DAY`, and it carries its author's own comment: not sure that
RADCON makes sense when the diurnal cycle and gravity are changed, so setting
them to use Earth values. An exoplanet-capable code, with a derived constant
that escaped the parameterisation, and the uncertainty left on the record rather
than resolved.

That is worth more than the finding itself, because it says the
inherited-constant class needs auditing regardless of what a code was built for.
So the obvious check was made here: at this planet's gravity and day length a
RADCON-like constant would be wrong by about 1.63x, which would be enormous.
**ExoPlaSim is clean.** `radmod.f90` computes the heating rate as
`-ga*(dflux(jlep)-dflux(jlev)) / (dsigma*dp*acpd*(1+ADV*dq))` and the layer
thickness as `gascon/ga`, using the model's own gravity and specific heat
throughout. A negative result, recorded because the check was cheap and the
failure mode is real enough that a purpose-built code fell into it.


## 10. cGENIE, first reading

*Read 2026-08-22 from `vendor/cgenie/` at the vendored commit. The tree was
vendored for OCN-3 and five rows wait on it; this is the reading that makes
their citations checkable, plus one thing none of them anticipated.*

### 10a. The citations verify, and two are stronger than the rows state

**OCN-19's dead duplication is real and slightly worse.**
`genie-goldstein/src/fortran/tstepo.F` declares six arrays at lines 32 to 35,
two of them full four-dimensional tracer arrays `ts_t1` and `ts_t2` plus
`ts1_t1`, `ts1_t2` and two density copies, and fills all six by copy at lines 62
to 68 on every ocean timestep. The routines that would consume them,
`tstepo_flux_t`, are COMMENTED OUT, as is the entire OpenMP block around them.
The only live call is `call tstepo_flux()`.

**OCN-19's "no disabled parallelism to re-enable" verifies and is stronger than
a judgement.** Both OpenMP blocks are commented out rather than merely inactive,
and the biogem one at `genie-biogem/src/fortran/biogem.f90:718-735` calls
`sub_wasteCPUcycles1` and `sub_wasteCPUcycles2`, which are real subroutines at
`biogem_lib.f90:2267` and `:2289`. The block was a scheduling experiment on
deliberately wasted cycles, not a parallelisation anyone disabled.

**OCN-20's ordering verifies at `invert.f:32`,** `k=i + j*n`, exactly as quoted.
What the row does not carry is the comment three lines above it:

    c NOTE: gfortran compiler will multiply flag as:
    c     'Array reference ... out of bounds' ... but all OK
    c     (I think ...)

An acknowledged, unresolved out-of-bounds in the routine OCN-20 proposes to
replace, hedged by its own author. That is a CORRECTNESS argument for the
replacement beside the complexity one, and this project has both the instrument
and the experience to settle it: `-ffpe-trap` is in the production flag line and
the checked profile carries `-fcheck=all` precisely because that class caught a
real out-of-bounds during the CLIM-40 fold work.

### 10b. cGENIE ships a PlaSim, and it is our PlaSim

`vendor/cgenie/genie-plasim/src/fortran` is a complete Planet Simulator:
`plasim.f90`, `plasimmod.f90`, `radmod.f90`, `oceanmod.f90`, `seamod.f90`,
`icemod.f90`, `landmod.f90`, `simba.f90`, `legmod.f90`, `fftmod.f90` and the
rest, identifying itself as PUMA 2.0, Version 16 Revision 4.

**Twenty-four of its twenty-five modules are also in
`vendor/exoplasim/exoplasim/plasim/src`.** ExoPlaSim adds to that set -- the
aerosol stack, glaciers, the carbon module, LSG, the MPI variants, the planet
configurations -- and removes nothing. The single cGENIE-exclusive file is
`geniemod.f90`.

That a PlaSim-to-cGENIE path exists is the PREMISE of OCN-17 rather than news:
that row exists because of the published `exoplasim_genie_regrid` work. What the
module comparison adds is how close the bundled one is -- written against the
same module names and the same `pumamod` state this project's model still has,
so it is a reference interface OCN-10 can read rather than a related project to
cite. It does not reopen OCN-3's offline decision, which rests on deep-ocean
spin-up cost and is untouched.

### 10c. The coupling contract, and a lead on OCN-17's multiplier

`geniemod.f90` is 71 lines declaring the exchange arrays and `fluxmod.f90` is
1,092 doing the work. What crosses, atmosphere to ocean: net solar, evaporation,
precipitation, runoff, and TWO PAIRS of wind stress. What comes back: sea
surface temperature, ice temperature, ice height, ice fraction and ice albedo.
Alongside: cell area, land-sea mask, an atmospheric pressure term, and CO2. One
array, `genie_dflux`, is annotated `!not used`.

Two details matter more than the list.

**The two stress pairs are a C-grid staggering, and cGENIE says so in its own
comments.** `genie-main/genie_loop_wrappers.f90` lines 34 and 36 read
`ocean_stressx2_ocn  surface wind stress (x) at u point` and
`ocean_stressx3_ocn  surface wind stress (x) at v point`. `genie.F` lines 342 to
360 accumulate both over `kocn_loop` steps and divide by the step count, so what
the ocean receives is a time mean of a stress evaluated at two staggered
locations.

That is a concrete mechanism for the knob OCN-17 exists to explain. The
published `exoplasim_genie_regrid` path instructs users to set `ea_11` and
`go_13` to 2.0 or 2.6 depending on ExoPlaSim version, with a matching change to
gas transfer. A regridding route that supplies one stress field where the ocean
expects two staggered ones, or that misses the `kocn_loop` averaging, would need
a compensating multiplier of about that size, and the fact that the scaling is
applied on BOTH the atmosphere and ocean side is consistent with an
interpolation mismatch rather than a physical correction. This is a hypothesis
with a mechanism and a place to look, not an explanation; OCN-17 still owns the
verdict.

**The exchange arrays are dimensioned 360.** `g_netsolar_plas`,
`g_evap_plas`, `g_precip_plas`, `g_runoff_plas` and all four stress arrays carry
a trailing dimension of 360, an Earth year of daily fields compiled in. Vesper's
orbit is 180.7 days. That is a calendar assumption in the coupling itself rather
than in a namelist, and OCN-10's contract has to state what replaces it.

### 10d. OCN-18's three questions, two answered and one needing a download

**The island bound is 5, not 10, and cGENIE's own configurations routinely
exceed it.** `genie-main/genie_control.f90:79` reads
`#define GOLDSTEINMAXISLES 5`. Counting the maximum island index in every
`.psiles` file shipped under `genie-paleo`, across 294 configurations:

| islands | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| configurations | 46 | 58 | 47 | 54 | 34 | 22 | 27 | 6 |

So 89 of 294, about thirty percent, need the macro raised above its default, and
none exceeds 8. Island count is therefore a per-geography property that drives a
COMPILE-TIME bound, and a new geography can require a rebuild before it will run
at all. That is a constraint on OCN-11's bathymetry contract as much as on
OCN-18: this world is about three quarters endorheic with a coastline nothing
has counted islands on, and the number is not knowable until a candidate
configuration is generated. OCN-20's instruction to cost the island machinery
first is sharpened by this, since `ratm`, `psisl`, `erisl` and `matinv_gold` all
scale with it.

**The resolution ceiling is declared and unexercised.** muffingen's own settings
files annotate `par_max_i` and `par_max_j` as `[1-72]`, four times the cell count
of the shipped grid. But of the shipped configurations 131 are 36 x 36, five are
18 x 18 and one is 12 x 12, and NOTHING ships above 36. So the range is a
declaration rather than a demonstration, which is the distinction OCN-18 exists
to make: the tier's attraction is its price and that price is a property of
36 x 36.

Two details worth carrying. `par_max_k=17` with `par_add_Dk=1`, so the shipped
depth structure is 16 levels plus one extra. And `par_A_frac_threshold=0.50`, a
land fractional area threshold identical in convention to this project's
`geography_land_threshold`.

**muffingen itself is not in this tree.** `genie-paleo` holds only its outputs,
its per-configuration settings files and its run logs; the generator is a
separate repository. Answering what it does above 36 x 36, and whether it runs
under Octave, needs that download. The settings and logs survive here only
because the earlier deletion took the PostScript plots and left everything else.

### 10e. GOLDSTEIN's non-dimensionalisation, which is half Earth's

`initialise_goldstein.F:370-398` sets the scales the whole ocean is
non-dimensionalised against. Some are configurable and some are not, and the
split is not where a reader would guess.

**Configurable, through `ini_gold_nml`:** `sodaylen` and `sidaylen`, the solar
and sidereal day lengths; `yearlen`; `nyear`; `par_dsc`, the depth scale; and
`par_dk`.

**Hardcoded in the source:**

    usc   = 0.05      ! velocity scale, m/s
    rsc   = 6.37e6    ! EARTH RADIUS
    gsc   = 9.81      ! EARTH GRAVITY
    rh0sc = 1e3       ! reference density
    cpsc  = 3981.1    ! seawater heat capacity

`rh0sc` and `cpsc` are properties of sea water and transfer. `usc` is a
modelling choice rather than a planetary constant. `rsc` and `gsc` are Earth's
and are not exposed anywhere.

**And four derived scales carry them into everything:**

    tsc     = rsc/usc                          ! every non-dimensional time
    rhosc   = rh0sc*fsc*usc*rsc/gsc/dsc        ! carries BOTH rsc and gsc
    opsisc  = dsc*usc*rsc*1e-6                 ! overturning streamfunction
    rfluxsc = rsc/(dsc*usc*rh0sc*cpsc)         ! heat flux scaling

At this planet's radius of 1.20 Earth and gravity of 1.306 Earth, `tsc`,
`opsisc` and `rfluxsc` are each wrong by a factor of 1.20 and `rhosc` by about
eight percent, in a model that is entirely non-dimensional. An overturning
reported through `opsisc` would be twenty percent low, and it would look
plausible.

**Rotation is the exception, and it fails open.** Somebody added sidereal-day
support in January 2024, and the code reads:

    if (abs(sodaylen-sidaylen).gt.0.001) then
       fsc = 4*pi/sidaylen
    else
       fsc = 2*7.2921e-5
    endif

The comment says the branch exists for backwards compatibility. What it means in
practice is that a configuration setting the solar and sidereal day lengths
equal -- an easy thing to do wrong, and exactly wrong for a rotating planet --
silently reverts the Coriolis scaling to EARTH's rate rather than failing. That
is the fail-open pattern this project's conventions forbid, in the one planetary
constant the source does parameterise.

There is also precedent in the file for a scale being wrong. Line 642 carries
`rma error: rhosc value in coeffs assumed dsc=4 km, corrected offline for IPCC
runs with GENIE-1`, so a scaling constant has already been wrong once here and
was patched outside the code.

None of this is an argument against cGENIE. It is the inventory OCN-12 exists to
take, arriving early, and the answer to OCN-3's question of how the candidate
consumes Vesper radius, gravity, rotation and calendar is now specific: calendar
and rotation yes, with a fail-open on rotation; radius and gravity no.

### 10f. muffingen, and the blank world

`references/muffingen/`, GPL-3.0, the grid generator whose outputs fill
`genie-paleo`. Three things bear on open rows.

**It does not need a GCM.** `EXAMPLE_BLANK.m` sets `par_gcm=''` with every
netCDF input name empty, so muffingen will build `.k1`, `.paths` and `.psiles`
from a supplied mask and topography alone. That is OCN-11's route: the Orogen
bathymetry rather than someone else's model output. `EXAMPLE_K2_supercontinent.m`
is the idealised-world precedent beside it.

**It writes the gas transfer parameter OCN-17 is asking about.** `muffingen.m`
emits `bg_par_gastransfer_a = 0.722` into the configurations it generates, and
the published regrid path instructs users to move that to 1.1 alongside the
wind-stress scaling. So the multiplier that row is chasing has a default here to
be compared against, which is a second place to look after the staggering in
section 10c.

**Octave is unresolved.** The tree is MATLAB `.m` with a LaTeX manual under
`DOCS/`, and a quick search finds no Octave statement either way. OCN-18 asks
whether it runs under Octave and that remains open; the manual is where to look.


## 11. PHYS-14 answered: the model ships the spectra and never re-weights them

*Computed 2026-08-22. The row was opened on an argument from the vegetation
precedent; this is the argument checked, and it is stronger than stated.*

### 11a. Snow reflects identically either side of the band split

`vendor/exoplasim/exoplasim/plasim/src/plasimmod.f90:428-432`:

    real :: dsnowalbmx(2) = 0.8
    real :: dsnowalbmn(2) = 0.4
    real :: dglacalbmn(2) = 0.6

Two-element band arrays initialised from a scalar, so element 1 and element 2
are equal. The simulation's snow reflects the same fraction in the visible as in
the near infrared, which is false on Earth and false here.

That is the SAME DEFECT the vegetation work already fixed.
`config/planet.yaml` records it in those words: codes 175 and 176 carry the
vegetation albedo "where one identical field used to assert that a canopy
reflects equally either side of the split". Vegetation was repaired and given a
band pair. Snow, glacier ice and sea ice were not.

It also means the k25v spectrum fix did NOTHING for snow. The two-band machinery
carries a star's spectral shape by moving flux between the bands; if both bands
hold the same number, moving flux between them changes nothing.

### 11b. The model ships the spectra its own constants came from

`vendor/exoplasim/exoplasim/surfacespecs.py` is a library of ECOSTRESS surface
reflectance spectra interpolated to the 965 wavelengths ExoPlaSim uses, which is
the grid `k25v.dat` is written on. It contains `iceblendmax`, `iceblendmin`,
`glacalbmin`, `seaicemax` and `seaicemin`, each documented as a re-weighted
snow, ice, water and frost blend chosen to yield a specific model reflectance.

Only `__init__.py` and `pRT.py` import it. Nothing weights those spectra against
the configured star to produce the model's albedo constants.

### 11c. The numbers, and a check that the method is right

Weighting each spectrum by a 5772 K Planck and by `k25v_hr.dat`, split at
0.75 um:

| spectrum | model key | constant | Sun b1 / b2 / broadband | K2.5V b1 / b2 / broadband |
| --- | --- | ---: | --- | --- |
| `iceblendmax` | `dsnowalbmx` | 0.8 | 0.984 / 0.644 / 0.816 | 0.983 / 0.585 / 0.736 |
| `iceblendmin` | `dsnowalbmn` | 0.4 | 0.509 / 0.305 / 0.408 | 0.508 / 0.273 / 0.362 |
| `glacalbmin` | `dglacalbmn` | 0.6 | 0.766 / 0.455 / 0.612 | 0.764 / 0.406 / 0.541 |
| `seaicemax` | `dicealbmx` | 0.7 | 0.894 / 0.530 / 0.714 | 0.892 / 0.473 / 0.631 |
| `seaicemin` | `dicealbmn` | 0.5 | 0.637 / 0.380 / 0.510 | 0.636 / 0.339 / 0.451 |
| `oceanblend` | `doceanalb` | 0.069 | 0.076 / 0.065 / 0.070 | 0.076 / 0.064 / 0.068 |

**The solar column reproduces every one of the six model constants**, the worst
by 0.016 and the ocean to three decimals. Six for six is the check that matters:
it establishes that these spectra ARE the provenance of the model's constants,
and it validates the integral before any conclusion is drawn from the K dwarf
column.

Read across and the effect splits cleanly in two. Band 1 barely moves, at most
0.002, because snow is flat and bright below 0.75 um. Band 2 falls by 0.04 to
0.06, because this star's flux sits deeper into the infrared WITHIN that band,
where snow darkens. On top of that the band-1 flux share falls from 0.556 to
0.406, moving weight onto the darker band. The first of those is what no fixed
per-band constant can capture; the second is what the two-band machinery would
capture if the bands differed.

**So the simulation's fresh snow is about 0.064 too bright, its old snow 0.038,
its glacier ice 0.059, its maximum sea ice 0.069 and its minimum sea ice
0.049**, and every one errs in the direction that weakens the ice-albedo
feedback at the cold end. Sea ice at 0.069 is the LARGEST single error in the
set, and it is the one surface with no configuration key at all.

**The ocean is the exception, and the exemption is worth stating.** `oceanblend`
gives 0.070 under the Sun and 0.068 under this star, because water's reflectance
is low and nearly flat across the split, 0.076 against 0.065. So the SPECTRAL
half of OCN-7's question is answered negatively: re-weighting the ocean albedo
for this host is worth 0.002 and does not earn its place. That row's
spatial-variation half is untouched.

Seven arrays carry this defect in total, and two of them, `dsnowalb` and
`doceanalb`, are commented "spectral weighted" in the source while holding one
value in both bands. The comment asserts what the value denies.

### 11d. What remains

The derivation is not the work any more; the numbers above are it, and
`iceblendmax` against `iceblendmin` is already the bracket PHYS-14 asked for
over grain size and impurity. What remains is the shape: a provenance-stamped
analysis product on the pattern of `analysis/vegetation_albedo.py`, an anchoring
convention for the band pair like `vegetation_albedo_bands`, keys in
`config/planet.yaml`, and a decision about the sea ice pair, which belongs to
`icemod` rather than `landmod` and has no config key at all today.


## 12. OCN-4's photon question, answered on both candidates

*Computed and read 2026-08-22 from `vendor/cgenie/genie-ecogem` and
`references/marbl`. OCN-4 asks for the K-star photon dependency of each
candidate; this is that, plus what the answer costs.*

### 12a. What this star is worth in photons

Integrating `k25v_hr.dat` against a 5772 K Planck over 0.2 to 5.0 um, with the
conventional marine PAR window of 400 to 700 nm:

| | Sun 5772 K | K2.5V | ratio |
| --- | ---: | ---: | ---: |
| fraction of flux in 400-700 nm | 0.369 | 0.313 | 0.849 |
| umol photons per J WITHIN the band | 4.566 | 4.673 | 1.023 |
| umol photons per W of total shortwave | 1.684 | 1.463 | **0.869** |

**Read the ratios and not the absolutes.** The absolute fractions depend on the
integration bounds and on a blackbody standing in for the Sun; the ratios are
robust to both.

The useful finding is the split. Almost the whole effect is the BAND FRACTION,
at 0.849, and the mean photon energy inside the band moves only 2.3 percent
because 400 to 700 nm is narrow enough that the spectral tilt within it barely
matters. So a single scalar PAR fraction is adequate in STRUCTURE for a
non-solar host. What is wrong is its value: this star delivers about 0.87 of the
Sun's photon flux per unit shortwave.

### 12b. The land has this and the ocean does not

`FRADPAR` is 0.4624, derived from the k25v spectrum over a 400 to 750 nm window,
and `biosphere/README.md` carries it against K2-18's 0.309. So the terrestrial
side of this question was answered when the spectrum was adopted.

Nothing has done it for the ocean, and the window is not obviously the same
number. The terrestrial window was widened to 750 nm on predicted K2V
photosystem peaks at 675, 711 and 746 nm, which REF-9 records as resting on four
references all marked `held`. A marine window has the opposite pressure on it:
water absorbs strongly beyond about 600 nm, so the light that reaches a
phytoplankton cell is bluer than the light at the surface, and a window widened
to 750 nm describes photons the water column has already removed. **Under a red
host the two effects compound rather than cancel** -- less blue arriving, and
what does arrive being weighted toward the part water takes first.

OCN-4 therefore cannot inherit `FRADPAR`. It has to decide the marine window on
its own evidence, and the number will not be 0.4624.

### 12c. How the two candidates take it, and MARBL is better here

**ECOGEM** computes `PAR(:,:) = PARfrac * dum_egbg_fxsw(:,:)` in `ecogem.f90`,
with `PARfrac` a scalar in `ini_ecogem_nml`. Configurable, which is the good
part. Attenuation is `k_w` for water and `k_chl` for chlorophyll, both scalars
in the same namelist, so the light field is SPECTRALLY FLAT below the surface as
well as at it. That is the same defect class as the surface albedo arrays in
section 11, arriving in the ocean: a single attenuation coefficient cannot
express that a red-shifted star's photons are removed faster.

That is OCN-6's question rather than a new one, and this identifies the specific
scalars it would replace. It also sharpens it: the error is one-signed, since a
solar-tuned `k_w` under a redder star OVERSTATES the euphotic depth.

**MARBL** takes PAR Column Fraction and shortwave as `interior_tendency_forcings`
entries, so the conversion is the driver's responsibility and never a constant
inside the model. For a non-solar host that is strictly better: the quantity
OCN-4 would have to override in ECOGEM is one MARBL asks for by design, and
`num_PAR_subcols` means it accepts a sub-column structure rather than one value.

On this axis MARBL wins, and it is the axis this project cares about most.
OCN-4's other criteria -- tracer and trait cost, the covarying trait space,
fixed Earth inventories -- are untouched by this and still favour ECOGEM's
trait-based ecology in principle. The comparison is now split rather than
settled.


## 13. The known-good suite, and an ordering constraint it imposes

*Read 2026-08-22. OCN-19 validates its threading against `genie-knowngood/`,
described there as per-component reference output for four configurations. It is
thinner than that, the comparator is better than expected, and the combination
creates an ordering nobody has stated.*

### 13a. What it actually is

Four netCDF files, one per directory:

| configuration | component | file | size |
| --- | --- | --- | ---: |
| `genie_eb_go_gs_ac_bg` | biogem | `fields_biogem_3d.nc` | 1.9 MB |
| `genie_eb_go_gs_el` | ents | `ents_yearav_0000000020.nc` | 64 KB |
| `genie_eb_go_gs` | goldstein | `gold_spn_av_0000000020_00.nc` | 312 KB |
| `genie_na_go_ni` | goldstein | `gold_spn_av_0000000010_00.nc` | 312 KB |

Three components, not four, since two configurations cover GOLDSTEIN. The
integrations are twenty and ten model years, so this is a short-run regression
surface rather than an equilibrated comparison -- which is adequate for what
OCN-19 wants it for, and worth knowing before it is asked to carry more.

**There is no ECOGEM reference.** If OCN-4 selects the trait-based tier, the
ecosystem component arrives with no known-good at all, and one would have to be
generated and blessed before any change to it could be validated. That is a cost
on ECOGEM's side of OCN-4's ledger that the row does not currently carry.

### 13b. The comparator is the right instrument and already exists

`genie-main/src/c/compare.cpp`, built as `nccompare.exe` and driven by
`compare-basic.sh`, takes `-r` as a relative tolerance IN ULPS and `-a` as an
absolute floor below which differences are ignored.

A tolerance in units in the last place is exactly the right measure for a
reordering or threading change, because it asks how far the result moved in
representable numbers rather than as a fraction. OCN-20 requires that the
tolerance and the acceptance comparison be fixed before an iterative scheme is
chosen; the instrument for that is in the tree and does not have to be written.

### 13c. The constraint: this suite validates EARTH, so optimise first

The GOLDSTEIN reference carries `opsi`, `opsi_a` and `opsi_p`, the global,
Atlantic and Pacific overturning streamfunctions. Those are reported through
`opsisc = dsc*usc*rsc*1e-6`, and section 10e establishes that `rsc` is Earth's
radius, hardcoded.

So every number in the known-good is a number at Earth's parameters. Change
`rsc` for this planet and the reference values move with it, in a scaling that
is not a simple factor across fields, and the suite stops being a comparison at
all.

**That orders the ocean work.** OCN-19's threading and OCN-20's solver
replacement are correctness-preserving changes at fixed physics, which is
precisely what a ULP-tolerance regression suite is for, and they can use it. The
radius and gravity work under OCN-3 and OCN-12 is not, and it destroys the
suite's applicability as its first act. Doing the optimisation first keeps a
validation surface that doing it second would not have.

If the order goes the other way, a replacement acceptance test has to be
declared before the parameters move, and it cannot be a comparison against
these files.


## 14. Which half of OCN-19 actually threads

*Static analysis 2026-08-22. OCN-19 names two threading targets, BIOGEM's
`do n=1,n_vocn` sweep and the 3-D loops in `tstepo_flux`, and treats them as one
piece of work. They are not comparable.*

### 14a. BIOGEM threads, on three checks

The loop in `genie-biogem/src/fortran/biogem.f90:754` is:

    do n=1,n_vocn
       call sub_box_remin_DOM(vocn(n),vbio_remin(n),loc_dtyr)
       call sub_box_remin_part(loc_dtyr,vocn(n),vphys_ocn(n),vbio_part(n),vbio_remin(n))
    end do

Every argument is indexed by `n`, `loc_dtyr` is a read-only scalar, and the
reduction back onto the grid, `bio_remin = bio_remin + fun_lib_conv_vocnTOocn(...)`,
happens AFTER the loop rather than inside it.

Both callees live in `biogem_box.f90`, at 161 and 1,208 lines. Checked for the
three things that would break it:

- no `SAVE` or `DATA` statements in either;
- no writes to module-scope arrays -- `bio_remin`, `bio_part`, `ocn`,
  `phys_ocn`, `bio_settle` are all reached through dummy arguments;
- **no initialised locals anywhere in `biogem_box.f90`**, which is the implicit
  SAVE that CLIM-51 catalogues, where a local declared with an initialiser is
  silently persistent and therefore silently SHARED under `-fopenmp`.

That last one is the striking result. CLIM-51 found 103 such declarations in
`plasim/src` and records that "gfortran reports nothing at any warning level
tried, so there is no flag that shortens this and it is reading". Across all of
BIOGEM and ECOGEM the same pattern appears 23 times, none of them in
`biogem_box.f90`:

| file | count |
| --- | ---: |
| `biogem_lib.f90` | 14 |
| `ecogem_lib.f90` | 4 |
| `ecogem.f90` | 2 |
| `biogem_data.f90`, `biogem_data_ascii.f90`, `ecogem_data.f90` | 1 each |
| **`biogem_box.f90`** | **0** |

Two honest caveats. That count is an UPPER BOUND: it does not separate
module-scope declarations, which are legitimately shared, from procedure-body
ones, which are the hazard, and CLIM-51's 103 counts only the second kind. And
this is one level into the call tree; the callees call further functions whose
files are among the 23.

### 14b. GOLDSTEIN does not, and the reason is structural

`tstepo_flux.F` includes `ocean.cmn` twice, and `ocean.cmn` declares
**44 common blocks**. GOLDSTEIN's entire state is global by construction:
`ocn_invars`, `ocn_vars`, `ocn_lego`, `ocn_islands` and forty more. The loops
themselves are Fortran-77 numbered nests over `j`, `i` and `l`.

Threading that is not the same job as threading BIOGEM. In BIOGEM the state
arrives as derived-type arrays indexed by the loop variable, so privacy is the
default and sharing is explicit. In GOLDSTEIN sharing is the default and every
variable the loop touches has to be proved read-only or made private, across
44 blocks.

### 14c. What this does to the row

OCN-19 should not treat its two targets as one decomposition. BIOGEM's sweep is
the tractable half and the static evidence says it threads; `tstepo_flux` is a
common-block privatisation exercise that happens to contain loops.

It also changes what the profile is FOR. OCN-19 orders the profile first so that
the half that pays can be identified, which is right. But if the profile says
`tstepo_flux` dominates, the answer is not "thread it" -- it is that the
threading route is expensive there and OCN-20's algorithmic route matters more.
The profile decides between two DIFFERENT kinds of work, not between two
candidates for the same kind.
