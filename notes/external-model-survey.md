# What other models in this class do, and what is worth taking

*Written 2026-08-22. Claims about this tree are from reading it; claims about
each external model are from its published description or its source, cited
inline. Costs quoted for other models are what their own authors state, on their
own hardware, and are not benchmarks run here.*

Vesper is a fictional planet and this project simulates it. What follows
compares this project's machinery against the Earth system models of
intermediate complexity, the exoplanet general circulation models, and two
recent projects that sit between them, so the comparison is not run again. It
reaches no adoption decision. What it produces is a list of formulations worth
taking, a list of candidates now closed, magnitudes for terms this project
carries as declared rather than measured, and corrections to several of this
project's own assumptions where an external implementation showed one to be
unsupported. Every finding is carried by a row in `TASKS.md`; the index below
is the check on that.

## Index

Each section and the task rows that carry it. A section with no row either
closes a candidate, records a negative result, or is framing; those say so.
This table is the check that nothing here is a finding without a home.

| # | what it settles | carried by |
| --- | --- | --- |
| 1 | The two families do not overlap, and this stack is in the gap | *framing; opens nothing* |
| 2 | Cost, and why the obvious comparison does not hold | `OCN-3` |
| 3 | PALADYN and its neighbours in CLIMBER-X | `CLIM-53`, `CLIM-63`, `DUST-14`, `GRAV-8`, `GRID-2`, `GW-26`, `LSHY-3`, `LSHY-5`, `LSHY-7`, `OCN-21`, `WET-10`, `WET-12` |
| 4 | BIG-MITgcm is the nearest peer, and two of its details are warnings | `CLIM-62`, `GRAV-6`, `OCN-11`, `PHYS-13` |
| 5 | CliMA is not a competitor and has one thing worth reading | `OCN-12` |
| 6 | Candidates this survey closes | *closes candidates; `notes/audits/dormant-exoplasim-modules.md` records them* |
| 7 | REF-10: the rest of the CLIMBER-X tree, triaged | *REF-10, now closed* |
| 8 | SPEEDY and SOCRATES, read for CLIM-61 | `CLIM-61` |
| 9 | The exoplanet models, read without a row to justify it | *parent of 15 and 21* |
| 10 | cGENIE, first reading | `OCN-10`, `OCN-11`, `OCN-12`, `OCN-17`, `OCN-18`, `OCN-19`, `OCN-20`, `OCN-3` |
| 11 | PHYS-14 answered: the model ships the spectra and never re-weights them | `OCN-7`, `PHYS-14` |
| 12 | OCN-4's photon question, answered on both candidates | `OCN-4`, `OCN-6` |
| 13 | The known-good suite, and an ordering constraint it imposes | `OCN-19`, `OCN-20`, `OCN-3` |
| 14 | Which half of OCN-19 actually threads | `OCN-19` |
| 15 | ExoRT is the wrong correlated-k for this planet, and the reason generalises | `CLIM-61` |
| 16 | OCN-4 compared: the split is scientific against engineering | `OCN-4` |
| 17 | Three published weathering schemes, shipped with their constants | `VOLC-9` |
| 18 | muffingen's wind, and a second mechanism for OCN-17's multiplier | `OCN-11`, `OCN-17` |
| 19 | What EMIC-class sea ice motion actually is | `OCN-21` |
| 20 | Six methane schemes and a burial bracket | `OCN-13`, `WET-10` |
| 21 | SOCRATES exposes the planet, not just the star | `CLIM-61` |
| 22 | An independent two-band ice albedo, and what it settles for PHYS-14 | `PHYS-14` |
| 23 | ClimaLand, and the land-column hypotheses LSHY registers | `GRAV-8`, `LSHY-1`, `LSHY-3` |
| 24 | The conversion FRADPAR left behind | `BIO-25` |
| 25 | The LPJ-GUESS Earth-constant scan, and where its results already live | *negative result: says where the scan already lives, so it is not repeated* |
| 26 | GEMlite: OCN-3's argument survives, and a fourth acceleration idiom | `CLIM-53`, `OCN-3` |
| 27 | cGENIE is flux-forced already, and OCN-5's requirement splits | `OCN-5` |
| 28 | Why OCN-5 has to be a loop, derived rather than asserted | `OCN-10`, `OCN-17`, `OCN-2`, `OCN-5` |
| 29 | Ocean albedo: three branches, and the two zenith forms are a free bracket | `OCN-7` |
| 30 | `genie-plasim`: the coupling contract, written down, and how it is afforded | `OCN-10`, `OCN-17`, `OCN-21`, `OCN-3` |
| 31 | VOLC-9's weathering columns, confirmed from use -- and the bracket is a different one | `VOLC-9` |
| 32 | Two one-bucket land surfaces that disagree by an order of magnitude | `LSHY-3`, `LSHY-4` |
| 33 | Egea's family is the shape vocabulary -- but it limits a different flux | `LSHY-3` |
| 34 | Anoxia: computed from air-filled porosity, or switched on by latitude | `BIO-29`, `WET-2`, `WET-6` |
| 35 | Two accelerators in one codebase, and only one of them is guarded | `CLIM-53`, `OCN-3`, `OCN-5` |
| 36 | A second TOPMODEL, and where its parameters come from | `GW-26` |
| 37 | Snow: a third density model, and an albedo predictor we do not have | `GRAV-8`, `PHYS-14` |
| 38 | The leaf-to-canopy bound is asserted, and it is not one-signed | `BIO-18` |
| 39 | Soil albedo does not know whether the soil is wet, and the blocker is the state variable | `DUST-17`, `PHYS-15` |
| 40 | OCN-12's inventory, geochemistry side: BIOGEM and SEDGEM | `OCN-1`, `OCN-12`, `OCN-3`, `OCN-4` |

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
point shallower than 20 m is set to 20 m. On a world whose terrain generator
produces closed basins in quantity -- the pre-carve build is 76 percent endorheic
by land drainage, which rule 9 makes a LIMIT rather than a state -- a preparation
step that silently removes enclosed water is not a detail. The carve loop exists
to open those basins, so the derived figure will be lower and is not yet known;
what matters here is that the class of feature exists at all and the step deletes
it irrespective of size. OCN-11 already requires bathymetric smoothing and
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
the filename as though it did would have been wrong.

CORRECTION: an earlier version of this section named `restore_salinity.f90`,
`flux_adj.f90` and `hosing.f90` as cGENIE's. They are CLIMBER-X's ocean files,
read in section 7a and misattributed here. What cGENIE actually has is section
27.

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

ExoRT ships FOUR band configurations, n28, n42, n68 and n84, chosen by
atmospheric regime, which looks at first like a cost ladder for CLIM-61.
Section 15 reads their gas coverage and withdraws that: they are aimed at
Archean and hydrogen-rich atmospheres, they carry no ozone and no N2O, and
`n68equiv` costs about three times `ga7`. The ladder is real and it is a ladder
for other planets than this one.

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
OCN-18: this world's coastline has never had islands counted on it, and the
number is not knowable until a candidate configuration is generated on an
accepted terrain -- which is after loop A closes, since carving changes the
coastline. OCN-20's instruction to cost the island machinery
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
for this host is worth 0.002 and does not earn its place.

**And it matters even less than that, for a reason that took a second reading of
`radmod` to find.** In the default configuration `doceanalb` never reaches open
water at all: `radmod.f90:2895-2904` overwrites the open-ocean term with a
zenith-dependent formula, leaving the constant in force only over the sea-ice
covered fraction. Section 29 has that in full, along with OCN-7's spatial half,
which it answers.

Seven arrays carry the flat-band declaration in total, and two of them,
`dsnowalb` and `doceanalb`, are commented "spectral weighted" in the source
while holding one value in both bands: the comment asserts what the value
denies. Of the seven it is the SNOW and ICE arrays that bind in practice --
`doceanalb` is a declaration the shortwave overrides, while `dsnowalb`,
`dsnowalbmx`, `dsnowalbmn`, `dicealbmx`, `dicealbmn` and `dglacalbmn` are all
used as written.

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

Section 40e puts a number on the second of those: only 0.382 of this star's flux
falls below 0.75 um against 0.537 of the Sun's, and water removes essentially
all the rest within half a metre. **Do not read that 0.712 as a rival to the
0.869 above.** They measure different things in different places. 0.869 is
PHOTON flux per unit shortwave inside 400-700 nm AT THE SURFACE; 0.712 is the
fraction of shortwave ENERGY still present BELOW the top metre or two. The first
sets what a surface-referenced PAR parameter is worth, the second how much
reaches anything living beneath the skin.

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


## 15. ExoRT is the wrong correlated-k for this planet, and the reason generalises

*Read 2026-08-22. Section 9a called ExoRT a second cost point for CLIM-61 with
structure SOCRATES lacks. Reading its gas coverage inverts that.*

### 15a. What it covers, and what it does not

`source/README` states the configurations plainly: `n68equiv`, the recommended
one, is CO2, H2O and CH4 from HITRAN2016 in 68 bins with equivalent-extinction
overlap; `n84equiv` is the same at 84; `n28archean` adds N2-N2, N2-H2 and H2-H2
CIA at 28 bins; `n42h2o` and `n68h2o` are N2 and H2O only.

The k-array roots in `n68equiv/radgrid.F90` confirm it: `kc` for line
absorption, then `kco2ch4`, `kco2co2`, `kco2h2`, `kh2h2`, `kh2ofrgn`,
`kh2oself`, `kmtckd`, `kn2h2`, `kn2n2`, `ko2co2`, `ko2n2`, `ko2o2`. A line
table plus a full collision-induced set.

**There is no ozone table and no N2O.** `kabs.F90` mentions neither.

That is disqualifying here rather than inconvenient. `atmosphere.ozone` is true,
`ozone_scale` is a derived 0.794, and `ozone_uv_weight` and
`ozone_visible_weight` are separately derived for this host, because the star's
ultraviolet is the reason `surface_uv_relative_to_earth` is 0.4. And
`pN2O_bar` is a measured 3e-07 with CLIM-42 having added a longwave band
specifically to carry CH4 and N2O. Adopting ExoRT would DISCARD both.

SOCRATES `ga7` lists water vapour, carbon dioxide, ozone, dinitrogen oxide,
methane, oxygen, sulphur dioxide and carbonyl sulphide -- every gas this world
has, and two spares.

### 15b. And it is more expensive, not less

`n68equiv` is 68 intervals at 8 Gauss points, so 544 monochromatic calculations
per column, and `n84equiv` is 672. SOCRATES `ga7` is 6 shortwave bands at up to
12 k-terms and 9 longwave bands, of order 170. The present scheme is 2 shortwave
bands and one longwave.

So the exoplanet-native code costs roughly three times the Earth code AND covers
fewer of this atmosphere's gases. On both axes that matter to CLIM-61 it loses.

### 15c. The reason, which is worth more than the verdict

ExoRT is built for the exoplanet cases the field actually studies: Archean
Earth, hydrogen-rich atmospheres, tidally locked worlds around M dwarfs. Those
need N2-H2 and H2-H2 CIA and have no ozone layer to speak of. Its band
configurations are named for exactly those regimes.

**Vesper is not that kind of exoplanet.** Its atmosphere is 1 bar of N2, O2 and
Ar with 450 ppm CO2, interactive water vapour and an ozone layer -- a composition
much closer to modern Earth than to any archetypal exoplanet target. What is not
Earth-like is the STAR.

That is the generalisation, and it should be applied when weighing any candidate
component: for this project, an Earth-built code with a substitutable stellar
spectrum will usually fit better than an exoplanet-built code aimed at exotic
compositions. Section 8b already showed SOCRATES treats the star as a data file.
This is the other half of the same argument, and it is why CLIM-61's candidate
order should be SOCRATES first rather than ExoRT.


## 16. OCN-4 compared: the split is scientific against engineering

*Read 2026-08-22 from `vendor/cgenie/genie-ecogem` and `references/marbl`.
Section 12 settled the photon axis; this is the rest of what OCN-4 asks for.*

### 16a. The trait space, where ECOGEM wins on the mechanism the row cares about

ECOGEM's plankton configuration is a data file. `8P8Z.eco` is a list of rows
reading type, size in micrometres, and a flag:

    Phytoplankton     0.60   1
    Phytoplankton     1.90   1
    ...  through 1900.00
    Zooplankton       0.60   1   ...

Eight phytoplankton on a log-spaced size spectrum from 0.6 to 1900 um, and eight
zooplankton. SIZE IS THE MASTER TRAIT and functional type the second axis, which
is exactly the declared trait space OCN-4 says trait-based ecology is preferred
for. And the cost is a file swap: the shipped set runs `1P` at one population,
`NPD`, `8P8Z` at sixteen, `8P7Z1F` and `8P7Z3F`, up to `32P32Z` at
sixty-four. No recompilation.

It also ships a MATCHED STRUCTURAL PAIR, `3Diat4ZP_PiEu.eco` against
`NoDiat4ZP_PiEu.eco`, which is a diatoms-on and diatoms-off comparison already
built. That is the shape this project's conventions ask for when pre-registering
a structural sensitivity, arriving free.

MARBL is the other thing. `autotroph_cnt` defaults to 1, is 3 under the CESM2
preset, and the settings file annotates it `cannot change : PFT_defaults ==
'CESM2'` -- so under the standard preset the count is LOCKED at three named
plant functional types. Composition is prescribed by PFT identity rather than
emerging from a trait space.

On the axis OCN-4 names as its reason for preferring trait-based ecology,
ECOGEM wins and it is not close.

### 16b. The engineering, where MARBL wins and it is also not close

MARBL ships `tests/regression_tests`, `tests/input_files/baselines`,
`tests/bld_tests`, and `MARBL_tools/code_consistency.py`. Its settings are YAML
with a JSON form, carrying per-configuration defaults, declared dependencies
such as `dependencies : base_bio_on`, longnames, units and datatypes, generated
by `MARBL_generate_settings_file.py`.

**And it ships a standalone driver**, `tests/driver_src` and `driver_exe`, so it
runs without an ocean model against prescribed forcing. For this project that is
worth more than it looks: OCN-3 names the aeolian pattern, downstream of one
climate and upstream of the next, as the only affordable coupling class, and a
component with its own offline driver is already shaped for it. OCN-4 could
exercise MARBL before any ocean host is selected.

ECOGEM has none of that. Section 13a records it has no entry in
`genie-knowngood`, and there is no separate test suite, no settings schema and
no offline driver. Its configuration is Fortran namelists and `.eco` text files.

Stoichiometry is a wash and worth recording as one, since OCN-4 asks what fixed
Earth relations survive. Both expose Redfield: MARBL as
`parm_Red_D_C_P`, `parm_Red_D_N_P` and six more in the YAML; BIOGEM as
`par_bio_red_POP_PON`, `par_bio_red_POP_POC` and others in `ini_biogem_nml`,
with `bio_part_red` a four-dimensional array so the ratios can vary spatially.
Neither buries them.

### 16c. What that leaves the row

The comparison does not resolve to a winner and should not be forced to. ECOGEM
has the ecology this project wants and none of the infrastructure; MARBL has the
infrastructure and a prescribed-PFT ecology.

Three things follow that the row can act on without deciding.

The photon axis in section 12 is not a tiebreaker either way now: `PARfrac` is a
namelist scalar in ECOGEM and settable, and MARBL's driver-supplied PAR is
cleaner but the difference is one number.

MARBL's standalone driver means it can be TRIED before an ocean host exists,
and ECOGEM cannot. That asymmetry is about sequencing rather than merit, and it
is the cheapest thing available: a scoping note could exercise one candidate for
real and only read the other.

And ECOGEM's cost ladder from one to sixty-four populations is the answer to
OCN-4's trait-count question, which the row currently frames as a single figure
of sixteen. Sixteen is the reference configuration, not the price.


## 17. Three published weathering schemes, shipped with their constants

*Read 2026-08-22 from `vendor/cgenie/genie-rokgem`.*

### 17a. What is there

rokgem's basic weathering is coarser than this project's: `par_weather_CaSiO3b`
and `par_weather_CaCO3` are GLOBAL rates in mol ALK per year, with
`par_k_Tb` and `par_k_Tg` as temperature exponents for basalt and granite and
`opt_weather_C_Si_bg` switching whether the two are distinguished at all. Against
a lithology map derived from Orogen's own generator, that is a step down.

Its spatially explicit mode is the interesting part. `par_weathopt` selects
`Global_avg`, `GKWM` or `GEM_CO2`, and `data/input` ships THREE published
parameterisations over lithology classes, each as a class list and a constants
table: Amiotte-Suchet 2003, Gibbs 1999, and GEM_CO2.

Their classes are carbonate, shale, sandstone, basalt, shield or granite, acid
volcanics, and ice.

### 17b. Where they agree and where they do not

Amiotte 2003 and GEM_CO2 share their first column exactly -- 1.586, 0.672,
0.152, 0.479, 0.095, then 0.272 against 0.222 for acid volcanics -- so they are
the same base rates. Gibbs 1999 is a different formulation on different units
entirely, 32000, 40000, 4000, 16000, 7900.

Where the two that share rates DIVERGE is the pair of columns that appear to
carry the carbonate against silicate split of CO2 consumption. For the carbonate
class GEM_CO2 reads 0.93 and 0.07; Amiotte reads 1.00 and 0.00. For sandstone
both read 0.125 and 0.875. The column meaning is inferred from position and the
schemes' published form rather than from a header, so it wants confirming before
being used.

**If that reading is right, the disagreement is exactly the thermostat.**
Silicate weathering is the long-term CO2 sink, because the carbon ends up in
marine carbonate; carbonate weathering returns its CO2 on precipitation and is
close to neutral over long timescales. So the fraction of weathering CO2
consumption attributed to silicate is what sets a carbonate-silicate
thermostat's strength, and two schemes sharing base rates still disagree on it.

### 17c. Why this project should care, and why it cannot simply adopt

`pedology/scripts/thermostat_efficiency.py` produces a thermostat number from
one parameterisation, and `config/planet.yaml` records that the outgassing
plausibility check behind the prescribed pCO2 came from
`weathering_fluxes.py`. There is no structural uncertainty attached to either,
and three published schemes disagreeing is a free estimate of what that
uncertainty is worth.

Adoption is a different matter and the class systems are why. rokgem works in
carbonate, shale, sandstone, basalt, shield and acid; this world's land is
playa clastic at about a quarter of it, continental clastic, schist, granite,
granodiorite, melange and evaporite. There is no land carbonate class here at
all -- LITH-26 places carbonate on the SHELF -- and playa clastic and evaporite
have no counterpart in any of the three schemes.

So the value is the BRACKET rather than the parameterisation: what does the
choice among published schemes cost the thermostat number, once each is mapped
onto the classes this world actually has. VOLC-9 owns that.


## 18. muffingen's wind, and a second mechanism for OCN-17's multiplier

*Read 2026-08-22 from `references/muffingen`.*

### 18a. What a from-scratch configuration actually produces

`EXAMPLE_BLANK.m` turns on seven generators: `opt_makemask`, `opt_maketopo`,
`opt_makeocean`, `opt_makerunoff`, `opt_makewind`, `opt_makealbedo` and the
sediment option. So muffingen does not just make the grid files; it makes a
runoff pattern, a wind field and an albedo file too.

**This project already has better versions of at least three of those.** Runoff
routing is `hydrography`'s, with drainage networks, catchments and a lake
solver. Albedo is the surface work of section 11 and `build_surface_albedo.py`.
Mask and bathymetry come from Orogen rather than being regridded from someone
else's model.

So OCN-11's contract is not "hand the export to muffingen and take the output".
It is: use muffingen for the grid-topology artifacts that only it knows how to
build, `.k1`, `.paths` and `.psiles`, and override the rest with what this
project already computes. Which of the seven are taken and which are overridden
is a decision the contract should state explicitly, because the defaults are all
ON.

### 18b. The wind it generates is synthetic and Earth-shaped

`source/make_grid_winds_zonal.m` describes itself as creating "a synthetic zonal
average windstress", and `par_tauopt` selects the regime: case 1 is low, annotated
modern Northern Hemisphere or paleo Eocene; case 2 is high, annotated water
world; case 3 detects a zonal gateway in the mask automatically. The file carries
separate zonal profile parameters for a world with land and for a water world.

For this planet that is the wrong object twice over. The profiles are fitted to
Earth's circulation, and this world's wind stress is something ExoPlaSim
computes at 30-hour rotation and 32 degree obliquity. `opt_makezonalwind` must
be off and the stress must come from the atmosphere.

### 18c. The gas transfer coefficient is not a second knob

`par_tauopt` does not only select the wind. In `muffingen.m` the same switch
writes `bg_par_gastransfer_a`, with the calibration chain in the comments: at
a = 0.310 the low case gives 0.0201 mol m-2 yr-1 uatm-1, the water world 0.0297
and the intermediate case 0.024903, leading to a = 0.722. All three branches
currently write 0.722, which is worth noting on its own since they read as
though they were meant to differ.

**So wind stress and gas transfer are coupled in the generator by design.** That
reframes OCN-17. The published regrid path instructs setting the wind-stress
scaling to 2.0 or 2.6 AND moving `bg_par_gastransfer_a` from 0.715 to 1.1, and
those are not two independent corrections. They are one adjustment about wind,
made in the two places muffingen itself couples.

It also says where the 0.715 came from: muffingen's own default is 0.722, so the
regrid path starts essentially at the generator's value and raises it.

That gives OCN-17 a second mechanism beside the C-grid staggering of section
10c, and this one is the more likely of the two. Replacing a synthetic zonal
stress with a real modelled one changes the stress magnitude and its
distribution at once, and the gas transfer coefficient is calibrated against
whichever it was. The question the row should ask first is what the ExoPlaSim
stress is being compared AGAINST -- muffingen's synthetic profile, or an
observational product -- because the multiplier is a ratio and the denominator
has not been identified.


## 19. What EMIC-class sea ice motion actually is

*Read 2026-08-22 from `vendor/cgenie/genie-goldsteinseaice`. OCN-21 was opened
today on the claim that sea ice dynamics needs wind stress and ocean currents.
That is half right, and the half that is wrong makes the eventual fix smaller.*

`tstepsic.f` announces itself at line 14 as a "2nd order explicit transport code
using **upper level ocean velocities**", and the code bears it out: `tstipsic.f`
builds its advective coefficients from `u(1,i,j)` and `u(2,i,j)`, the ocean's
top-level velocity components, with `diffsic` as a diffusion term beside them.
Ice area and thickness are carried as `varice` and moved by those coefficients,
with `par_sica_thresh` and `par_sich_thresh` gating the doorway cases.

Notably, a separate ice velocity field `uice` exists and is COMMENTED OUT at
every use site. The live code advects ice with the water underneath it.

**So this is not dynamic sea ice in the rheological sense.** There is no ice
momentum equation, no internal ice stress, no wind stress applied to the floe.
It is advection by ocean surface currents plus a diffusion term.

That is a ladder rather than a binary, and it is worth stating because the two
ends cost very different amounts:

| | what it is |
| --- | --- |
| this project | thermodynamic only; ice forms and melts in place |
| cGENIE | advection by upper-ocean velocity plus diffusion |
| CLIMBER-X | `sic_dyn.f90` and `transport_sic.f90`, a full dynamic-thermodynamic scheme |

**OCN-21's block is confirmed and its shape corrected.** The requirement is
narrower than wind stress and currents: cGENIE needs only the ocean's surface
velocity, and this world's ocean is a slab with no velocity field at all, so
there is literally nothing to advect with. `u(1,i,j)` has no counterpart here.

But the fix, when the block lifts, is smaller than "implement dynamic sea ice"
suggests. It is an advection term on two prognostic fields using a velocity the
ocean would already be producing. That belongs in OCN-21 because it changes what
the row is estimating, and it strengthens the case for declaring the term in the
error budget now: the eventual correction is cheap enough that the declaration
is not deferring an expensive decision.


## 20. Six methane schemes and a burial bracket

*Read 2026-08-22 from `vendor/cgenie/genie-atchem` and `genie-sedgem`.*

### 20a. WET-10's lifetime half has a route, not just a floor

`atchem_box.f90` carries SIX methane oxidation schemes:
`sub_calc_oxidize_CH4_default`, `_schmidt03`, `_claire06`, `_claire06_fixed`,
`_claire06H` and `_goldblatt06`, plus `sub_calc_reduce_N2O_lifetime` beside
them.

The default is Osborn and Wigley (1994), a power law in methane concentration:

    loc_tau = const_pCH4_oxidation_tau0 * (loc_CH4/const_pCH4_oxidation_C0)**const_pCH4_oxidation_N

with the comment "omitting [OH], [NOx] etc etc". The OH chemistry is folded into
three constants, and the concentration is clamped at half the reference value
because that is where the calibration curve ends -- Osborn and Wigley's curve
runs to four times it at the other end.

**`schmidt03` is the same functional form with the constants EXPOSED.** It reads
`par_pCH4_oxidation_tau0`, `C0` and `N` from the namelist where the default uses
hardcoded `const_` values, and its header names its calibration: a fit to the
2-D photochemistry model of Schmidt and Shindell (2003).

That is the route for WET-10 rather than merely a floor. The row needs an
oxidant and lifetime calculation for this star; `schmidt03` needs three numbers
from a photochemistry calculation, and `config/planet.yaml` already records that
Rugheimer et al. (2013) computed the photochemical steady state at this
effective temperature -- the same class of calculation Schmidt and Shindell's
model does. So the shape is: take `schmidt03`'s form, re-derive its three
constants for this host.

All the schemes conserve explicitly, which is the other thing WET-10 wants:

    O2  -= 2.0 * fracdecay * CH4
    CO2 += fracdecay * CH4          (with 13C and 14C carried through)
    CH4 *= (1.0 - fracdecay)

Methane oxidation consuming O2 two-for-one, the carbon arriving as CO2, and the
isotopes following. And it operates per grid cell rather than as a global box.

The `claire06` and `goldblatt06` families are the exception and should be left
alone: they take a whole-atmosphere conversion factor with no grid indices, and
they are Archean anoxic-atmosphere formulations. That is section 15c's lesson
again -- built for a kind of planet this one is not.

### 20b. Burial has published alternatives, like weathering did

`sedgem_lib.f90` selects diagenesis by string: `par_sed_diagen_CaCO3opt`,
`par_sed_diagen_opalopt` and `par_sed_diagen_Corgopt`, and the source ships four
sediment flux implementations beside them --
`sedgem_box_archer1991_sedflx.f90`, `sedgem_box_ridgwell2001_sedflx.f90`,
`sedgem_box_ridgwelletal2003_sedflx.f90` and `sedgem_box_benthic.f90`.

Organic carbon burial is fractional and redox-dependent:
`par_sed_diagen_fracCpres_ox`, `_anox` and `_eux` give separate preservation
fractions under oxic, anoxic and euxinic conditions, with
`par_sed_diagen_fracCpres_scale` for a Dunne scheme, and
`par_sed_huelse2017_sim_P_loss_pres_fracC` ties phosphorus loss to organic
carbon burial.

OCN-13 requires burial in its marine ledger and phosphorus among its species.
This is the same shape as section 17's weathering finding: several published
parameterisations shipped with their constants, so the bracket over scheme
choice is available rather than needing to be constructed. The redox dependence
matters here specifically, because OCN-13 also asks for anoxia and a redox
ledger, and these three fractions are where those two requirements meet.


## 21. SOCRATES exposes the planet, not just the star

*Read 2026-08-22 from `references/socrates/src/modules_core/rad_ccf.F90`.*

Section 8b established that SOCRATES treats the STAR as a data file. It treats
the PLANET the same way, and this is the only external code read today that gets
that structurally right.

Six constants are declared `PROTECTED` rather than `PARAMETER`, so they are
read-only to users of the module and writable within it:

    planet_radius   = 6.37E+06
    grav_acc        = 9.80665
    mol_weight_air  = 28.966e-03
    r_gas_dry       = 287.026
    cp_air_dry      = 1.005e+03
    repsilon        = 18.0153 / 28.966

And the module contains `set_socrates_constants(iu_nml)`, whose own comment
reads "Read planet-specific constants from namelist":

    NAMELIST /socrates_constants/ &
      planet_radius, mol_weight_air, grav_acc, r_gas_dry, cp_air_dry
    READ(NML=socrates_constants, UNIT=iu_nml)
    repsilon = molar_weight(ip_h2o)*1.0E-03_RealK / mol_weight_air

**`repsilon` is derived after the read rather than exposed**, so the water to air
molecular weight ratio cannot drift out of consistency with the air molecular
weight somebody set. That is the discipline `config/planet.yaml` already uses,
where mass follows from gravity and `derive()` raises if the two disagree,
arrived at independently.

Against everything else read today the contrast is sharp. CLIMBER-X puts Earth's
radius, gravity and even the Sun's visible share in a compile-time `parameter`
block, section 7a. GOLDSTEIN hardcodes `rsc` and `gsc` inside a subroutine and
derives four scales from them, section 10e. Isca uses a runtime namelist for
most of it and still let `RADCON` through as a `parameter`, section 9d. SOCRATES
is the one that exposes the set, protects it from outside writes, and derives
the dependent member.

For this planet only two of the five actually move. Radius goes to 7.645e6 and
gravity to 12.81; `mol_weight_air`, `r_gas_dry` and `cp_air_dry` are essentially
Earth's already, because this atmosphere is 78 percent N2, 21 percent O2 and
0.93 percent Ar, which is Earth's composition to three figures.

**So CLIM-61's candidate is now well-supported on five separate axes**, none of
which needed a run to establish: it covers every gas this world has, section
15a; it costs about an order of magnitude over the present scheme and roughly a
third of ExoRT, sections 8b and 15b; the stellar spectrum is a BT-Settl file and
the reweight is 46 lines, section 8b; the host-model interface is bounded at
3,452 lines by Isca's existing integration, section 9c; and the planetary
constants are a namelist with derived consistency. What remains genuinely open
is the cost per timestep in this model, which is the measurement the row is
ordered after the SHTns work to make.


## 22. An independent two-band ice albedo, and what it settles for PHYS-14

*Read 2026-08-22 from
`references/exocam/cesm1.2.1/configs/cam_land_fv/SourceMods/src.clm/SurfaceAlbedoMod.F90`.*

ExoCAM modified CLM's land ice albedo, and the modification is instructive twice
over. The header reads:

    ! The CLM default albice values are too high.
    ! Full-spectral albedo for land ice is ~0.5 (Paterson, Physics of Glaciers, 1994, p. 59)
    ! This is the value used in CAM3 by Pritchard et al., GRL, 35, 2008.

      real(r8), public  :: albice(numrad) = &   ! albedo land ice by waveband (1=vis, 2=nir)
                           (/ 0.80_r8, 0.55_r8 /)

**First: the two bands hold different numbers.** An exoplanet GCM carrying a
two-band surface albedo uses 0.80 in the visible and 0.55 in the near infrared,
not one value twice. That is independent confirmation that the seven flat arrays
in `plasimmod.f90` are anomalous rather than conventional, which is the claim
failure mode class 31 rests on.

**Second, and more useful: it shows the anchoring convention.** The values are
chosen against a stated FULL-SPECTRAL target of about 0.5 from Paterson, with a
precedent named. That is the same convention `vegetation_albedo_bands` already
uses here, where the band pair is "anchored so their flux-weighted combination
is the broadband value".

So PHYS-14's remaining shape question is settled by precedent on both sides:
this project's own vegetation work and an independent implementation reached the
same rule. The band pair is anchored to reproduce the broadband value, and the
broadband value is the one derived from the spectrum.

The numbers are close enough to be a sanity check and different enough to be
worth reading carefully. Section 11's `glacalbmin` derivation gives 0.766 and
0.455 under the Sun, against CLM's 0.80 and 0.55; the visible agrees well and
the near infrared sits lower. That is expected rather than troubling --
`glacalbmin` is ExoPlaSim's GLACIAL MINIMUM blend, meaning old or dirty ice,
while Paterson's 0.5 is a general glacier figure. The two are not the same
surface, and the ordering is the right way round.

Under `k25v` the same blend gives 0.764 and 0.406. The visible barely moves and
the near infrared falls by 0.049, which is section 11's finding arriving from a
second spectrum library.


## 23. ClimaLand, and the land-column hypotheses LSHY registers

*Read 2026-08-22 from `references/climaland/`, Apache-2.0, pulled because the
biosphere and land half of the task list -- BIO, ANUT, BVOC, DEMO, FIRE, PCAR,
PLHY, SDEC, EFOR, some ninety open rows -- had no external comparison at all.
Every other tree in this survey is a climate or ocean model.*

### 23a. It ships LSHY-3's hypotheses as interchangeable components

`src/standalone/` contains `Bucket`, `Soil`, `Snow`, `Vegetation`,
`SurfaceWater` and `InlandWater` as peer standalone models. LSHY-3 asks to
"register the current bucket, a parsimonious multilayer scheme and
gradient-driven flow as explicit hypotheses". ClimaLand ships the first and the
third of those as swappable components against a shared interface, which is the
SUMMA pattern of Clark et al. (2015) implemented -- and Clark is already in
`references/INDEX.md`'s LSHY table as the argument FOR registering hypotheses
this way.

`Soil/rre.jl` is the Richards equation. `Soil/energy_hydrology.jl` is the
phase-aware thermal and hydraulic column LSHY-5 describes. `Soil/Runoff` is
LSHY-2's partition. `Vegetation/canopy_energy.jl` and
`canopy_turbulent_fluxes.jl` are LSHY-4's fast vegetation-climate water loop.
`Soil/Biogeochemistry` is SDEC's ground. `Vegetation/pfts.jl` is PCAR-5's trait
registry.

### 23b. The retention closure is a choice this project has not recorded as one

`Soil/retention_models.jl` declares an abstract type with two implementations:

    export AbstractSoilHydrologyClosure, vanGenuchten, BrooksCorey

and `RichardsParameters` is parameterised over the closure type, so the choice
is made at construction rather than compiled in.

**That names a decision LSHY-1 currently leaves implicit.** LSHY-1 asks for a
retention curve among the properties in its contract and says to compare
pedology's AWC, ExoPlaSim's `dwmax` and LPJ-GUESS's Cosby derivation. Cosby
(1984) is a Clapp-Hornberger form, and PALADYN uses Clapp-Hornberger too, per
section 3. So this project sits entirely in one of the two families without
having written down that a family was chosen.

The two are not interchangeable in their tails. Brooks-Corey has a
discontinuous air-entry point and a power-law tail; van Genuchten is smooth
through saturation. Which matters most exactly where this world is unusual: a
quarter of its land is playa clastic in the pre-carve build, so wetting and
drying through near-saturation is a common state rather than an edge case.

That qualifier is load-bearing and I first wrote this without it. Playa extent is
NOT independent of endorheism: Orogen assigns basin fill on
`when: (c) => !!c.endorheic`, so opening a basin removes the cover and the
fraction moves with the carve list exactly as the drainage figures do. LITH-27
owns rechecking it. The argument survives in a weaker form -- this terrain
produces near-saturation surfaces as a class -- and not as a claim about a
settled area.

That is not an argument for switching. It is an argument that LSHY-1's contract
should carry the closure as a NAMED property with its family stated, so that a
consumer knows which one it is reading and a later comparison has something to
vary.

### 23c. What this does not supply

ClimaLand is CMIP-class and GPU-oriented, and section 2 established that this
project cannot use models in that throughput band. Nothing here is a candidate
for adoption; `Bucket` and `Soil` are read for their STRUCTURE, which is what
LSHY-3 needs, and PALADYN remains the EMIC-class implementation reference
because it runs at a comparable cost.


### 23d. Snow density, and a moderation of GRAV-8

GRAV-8 was opened on PALADYN's self-loading compaction term being linear in `g`
while ExoPlaSim carries a constant. A third implementation changes how that
should be read.

ClimaLand declares `AbstractDensityModel` in `Snow.jl` and ships exactly one
implementation, `MinimumDensityModel`, which computes

    ρ_snow = ρ_min * (1 - q_l) + ρ_l * q_l

a minimum dry-snow density blended toward liquid water density by the liquid
fraction. There is no self-loading compaction and gravity does not appear.

So across three land models the treatments are:

| | snow density |
| --- | --- |
| ExoPlaSim | constant 330 kg/m3 |
| ClimaLand | minimum density blended by liquid water content, no compaction |
| PALADYN | prognostic: Anderson fresh-snow temperature dependence plus Kojima self-loading, LINEAR IN g |

**Two of the three omit compaction, so ExoPlaSim's constant is conventional at
this complexity level rather than anomalous.** That is a real difference from
section 11's flat band arrays, where the structure existed and the value denied
it; here the structure is absent everywhere and only PALADYN builds it.

It does not make the term unimportant here, and the reason is the same one that
raised the question. A term that can be omitted on Earth is one whose
Earth-magnitude somebody judged small. That judgement does not transfer
automatically at 1.306 times Earth's gravity, and nothing in these three models
was asked to.

So GRAV-8's question is not "why is this missing", since it is missing almost
everywhere. It is whether the Earth-calibrated judgement that it is negligible
survives this planet's gravity -- which is what the row already asks for, a
price rather than a fix.

One detail worth carrying into that pricing: no model does all of it. ClimaLand
captures wet-snow density through liquid content, which PALADYN explicitly
neglects; PALADYN captures fresh-snow temperature dependence and self-loading,
which ClimaLand omits. The union of the two is what a complete treatment would
be, and nobody at this complexity level has built it.


## 24. The conversion FRADPAR left behind

*Computed 2026-08-22. Found by asking what ClimaLand's `λ_γ_PAR` corresponds to
here, which turned out to be a constant nobody had moved.*

ClimaLand converts absorbed PAR energy to mol photons with a representative
wavelength, `λ_γ_PAR`, set per experiment and defaulting to 500 nm. LPJ-GUESS
has the same quantity as `CQ`, and `vendor/lpj-guess/modules/canexch.h:42` reads:

    /// conversion factor for solar radiation at 550 nm from J/m2 to mol_quanta/m2
    const double CQ = 4.6e-6;

`canexch.cpp` uses it in every photosynthesis expression beside `apar`, and
`driver.cpp` builds `apar` as `rad * FRADPAR`. So `FRADPAR` and `CQ` convert
irradiance into photon-limited assimilation TOGETHER.

**`FRADPAR` is derived for this star. `CQ` is a `const` in a header, absent from
`vesper.h`, untouched by `build_vesper_header.py`.**

Integrating `k25v_hr.dat` against a 5772 K Planck:

| window | Sun | K2.5V |
| --- | ---: | ---: |
| 400-700 nm | 4.567e-6 | 4.673e-6 |
| 400-750 nm | 4.743e-6 | 4.864e-6 |

The solar 400-700 value reproduces the shipped 4.6e-6, which is the check that
the integral is the right one before anything is concluded from the rest.

Two errors compound and the larger is this project's own.
`biosphere/notes/productivity-prediction.md` records the photosynthetic window
being widened from 400-700 to 400-750 nm on Lehmer et al. (2021)'s predicted K2V
peaks at 675, 711 and 746 nm, moving `FRADPAR` from 0.3963 to 0.4624, a rise of
16.8 percent. `CQ` stayed at a value defined for the narrower and bluer band, so
the joules that widening ADDED are converted at a wavelength belonging to a
window that no longer applies.

Against FRADPAR's actual window and this star, `CQ` should be 4.864e-6 rather
than 4.6e-6 -- about 5.7 percent low, of which roughly 3.1 points is the window
mismatch and 2.6 the stellar shift. One-signed: more photons per joule than
assumed, so absorbed photon flux and therefore assimilation are understated.

**BIO-25 owns this.** `biosphere/notes/implicit-earth-assumptions.md` finding 5
is headed "The PAR correction changes energy but not photons per joule" and
specifies the fix: generate a spectrum- and window-weighted `VESPER_CQ` beside
`FRADPAR` and require the combined photon supply to reproduce the registered
productivity calculation. The plant-physiology audit repeats it. PCAR-11 was
opened against the same mechanism and withdrawn as redundant.

What was NOT already recorded is the magnitude, and that is folded into BIO-25:
about 5.7 percent low, decomposing into roughly 3.1 points of window mismatch
and 2.6 of stellar shift. The window half would exist on Earth's spectrum too,
which matters because it means the pair has to be recomputed over a COMMON
window rather than merely re-starred.


## 25. The LPJ-GUESS Earth-constant scan, and where its results already live

*2026-08-22. Kept so the scan is not repeated.*

Scanning `vendor/lpj-guess` for constants with an Earth basis, in the way
section 24 found `CQ`, turned up `PEATLAND_WETLAND_LATITUDE_LIMIT = 40.0` and
nothing else the project did not hold. **Do not re-run this scan expecting an
inventory.** `biosphere/notes/wetlands-peat-methane-audit.md` already has it,
under the heading "Two unrelated models selected by Earth latitude": it names
`is_highlatitude_peatland_stand()` as `PEATLAND && lat >= 40.0` and
`is_true_wetland_stand()` as `PEATLAND && lat < 40.0`, states the signed-latitude
consequence outright -- "a Vesper cell at 60 S follows the low-latitude
inundated-soil path" -- and records that the branch "controls hydrology,
vegetation stress, decomposition and methane, not merely a default parameter".
It prompted `implicit-earth-assumptions.md` finding 8 and BIO-29; WET-3 and
WET-6 own the consequences.

What the inventory did NOT carry was a magnitude, and section 34 supplies one:
on the slow and passive soil carbon pools the two branches differ by 0.36
against 0.025, a factor of 14.4, which is the pair that decides whether peat
accumulates at all.

**One method note that generalises**, because it cost a wrong claim here.
`is_true_wetland_stand()` reads as a property of the stand and is a property of
the stand AND its latitude. **A call site cannot be audited for latitude
dependence by looking at the call site**; the predicate has to be opened. That
is this project's convention about checking claims against the artifact rather
than the documentation, applied to a PREDICATE rather than a number, and BIO-29's
inventory has to be run that way to be worth anything.


## 26. GEMlite: OCN-3's argument survives, and a fourth acceleration idiom

*Read 2026-08-22 from `vendor/cgenie/genie-gemlite` and `genie-main/genie.F`.
Probed because OCN-3's choice of offline coupling rests entirely on deep-ocean
spin-up cost, and a module named `gemlite` beside `goldlite` and `ocnlite`
looked like it might undercut that.*

**It does not, and the reason is worth stating precisely.** `genie.F` cycles it:
`if (mod(koverall, kgemlite*kocn_loop) .eq. 1)`, running `gem_notyr` ordinary
years and then switching to a GEM cycle. On entry the cycle calls
`gemlite_cycleinit_wrapper`, then `cpl_comp_gemglt_wrapper` to "copy current
state of tracer arrays to GEMlite", then `gemlite_climate_wrapper` to "copy
climate state variables (here: sea-ice)".

So the CLIMATE STATE IS COPIED AND HELD while the geochemistry advances. GEMlite
is offline-tracer acceleration: it gets carbonate chemistry to equilibrium
cheaply once the circulation is settled. It does nothing for the circulation
itself, which is what OCN-3 means by equilibrating a deep ocean taking thousands
of model years. That row's argument for offline coupling is unaffected.

Its own header records the limit of the trick: surface fCO2 is rapidly
equilibrated with the atmosphere, so a pulse of CO2 emitted to the atmosphere is
NOT handled correctly, because GEMlite needs disequilibrium between ocean
surface and atmosphere to work. An accelerator that is wrong for transients.

### 26a. Four accelerators, three idioms

This tree now has four devices of the same family in reach, and they are not the
same shape:

| device | idiom | what it holds |
| --- | --- | --- |
| `NCVEG` | internal timestep multiplier | vegetation carbon advanced against fixed climate |
| PALADYN equilibrium spinup | internal timestep of 1000 years, possible because the components are fully implicit | annual cumulated NPP and litterfall |
| `newsnow` | extrapolate a measured 5-year tendency forward by a declared interval | snow where it is shrinking or has persisted |
| GEMlite | CYCLE: N ordinary years, then M accelerated years, repeat | climate state copied in and held |

GEMlite's is the most conservative of the three idioms, because it returns to
the full model periodically rather than running the accelerated component to
convergence in one stretch. That is a design worth knowing about for CLIM-53,
which asks not whether an accelerator is correct but what it is VALID FOR: a
cycling accelerator is valid under weaker assumptions than a single long jump,
since the full model re-establishes the state it was drifting from.

And all four state their invalidity somewhere, which is the property CLIM-53 is
really asking to reproduce. PALADYN's cannot be applied to processes
intrinsically out of equilibrium. GEMlite's cannot handle a transient
disequilibrium pulse. `newsnow` only extrapolates where snow already exists and
is shrinking or has persisted a year. The device is never the deliverable; the
statement of what it is valid for is.


## 27. cGENIE is flux-forced already, and OCN-5's requirement splits

*Read 2026-08-22 from `vendor/cgenie/genie-goldstein`. OCN-5 requires flux
forcing with no temperature restoring, so the question is what the candidate
does by default.*

**There is no sea-surface temperature or salinity restoring.**
`surf_ocn_sic.F` computes surface FLUX terms -- `fx0o` for heat into the ocean,
`evap`, `fxsen` and the sea-ice counterparts -- under a comment reading "main
i,j loop to compute surface flux terms". The ocean takes fluxes from whichever
atmosphere it is coupled to.

The one name that looks like restoring is not. `rel` sits in
`common /ocn_relax/` and `velc.f:103` uses it as

    u(1,i,j,k) = rel*u1(1,i,j,k) + (1.0 - rel)*u(1,i,j,k)

which is a numerical under-relaxation blending the new velocity with the
previous iterate. A solver control, not a physical term, and it should be
recorded as such so nobody retires it as a restoring.

The genuine explicit control is `get_hosing.F`, an additional freshwater
forcing with `hosing`, `hosing_trend` and `nyears_hosing` in `ini_gold_nml`. It
is called unconditionally from `goldstein.F:259` and its output is added to the
freshwater flux at line 297, so it is always in the path and inert only by
parameter value. That is exactly the class OCN-5 wants named rather than
removed, and it needs stating in the provenance chain because a zero default is
not the same as an absent term.

### 27a. The requirement splits, and only half of it is about the ocean model

OCN-5's warning is that "an ocean model forced by a climatology a zero-transport
slab produced, and restored toward that climatology's own sea-surface
temperature, returns a transport of zero by construction". That is two failures
joined by an "and", and cGENIE's construction only rules out the second.

There is no restoring, so the ocean cannot be pinned to a prescribed SST. But
nothing about the model prevents the first: if the heat and freshwater fluxes
handed to it are derived from a climatology a zero-transport slab produced, the
ocean will reproduce the transport that climatology implies, which is zero. No
property of the candidate protects against that, because the defect is in what
this project supplies rather than in what the model does with it.

So the row's requirement should be read as two, with different owners. "No
temperature restoring" is satisfied by the candidate and can be checked once.
"Flux forcing that does not encode the answer" is this project's, is not
checkable by inspecting cGENIE, and is the one that needs the declared
procedure the row asks for.


## 28. Why OCN-5 has to be a loop, derived rather than asserted

*Read 2026-08-22 from `vendor/cgenie/genie-embm`. Section 27 left OCN-5's real
risk as "flux forcing that encodes the answer" without saying what the mechanism
is. EMBM says it.*

### 28a. cGENIE tunes its atmosphere to make room for its ocean

EMBM moves heat and moisture by advection plus diffusion. `tstipa.f` builds the
coefficients as `betaz(l)*uatm(1,i,j)` zonally and `betam(l)*uatm(2,i,j)`
meridionally, with `diffa(l,1,j)` and `diffa(l,2,j)` beside them, indexed by
tracer and latitude. And `ini_embm_nml` exposes the lot: `diffamp(2)`,
`diffwid`, `difflin`, `betaz(2)`, `betam(2)`, with the latitudinal shape set by
`diffend = exp(-(0.5*pi/diffwid)**2)`.

So the atmosphere's poleward transport is a PARAMETERISED object with an
amplitude, a latitudinal width, a linear term and per-tracer advection scalings.
In a coupled cGENIE run those are calibrated so that EMBM and GOLDSTEIN TOGETHER
produce a sensible total transport. The partition between them is a tuned
quantity.

### 28b. This project cannot do that, and that is the whole problem

Poleward heat transport is shared between atmosphere and ocean. In an ExoPlaSim
run on a slab with `nhdiff` off, the ocean transports nothing, so the ATMOSPHERE
carries the entire load -- and it does so through resolved primitive-equation
dynamics, not through a coefficient anybody can turn down.

Hand that run's surface fluxes to a dynamic ocean and the ocean transports heat
too. The total is then an atmosphere already carrying the full load plus an
ocean carrying more. There is no `betaz` to reduce, because ExoPlaSim's
transport is emergent.

That is the mechanism behind OCN-5's warning, and it is sharper than the row's
own statement. The row says a restored ocean returns zero transport by
construction. The flux-forced failure is the opposite sign: an ocean forced by a
slab-derived climatology is being handed a flux field computed under the
assumption that it does nothing, and acting on it OVER-transports.

### 28c. But the published path does not do that, and the difference is an
### architectural fork nobody has stated

Section 28b cannot be the whole story, because the published
ExoPlaSim-to-cGENIE work exists and functions. The resolution is that it uses a
DIFFERENT ARCHITECTURE from the one OCN-10 is specifying.

EMBM is still running in it. `ea_` prefixed parameters appear 13,688 times
across the shipped configurations against 7,314 `go_` and 2,759 `gs_`, and the
regrid README instructs setting `ea_11`, an EMBM parameter. And EMBM's namelist
reads its winds from files: `xu_wstress`, `yu_wstress`, `xv_wstress`,
`yv_wstress`, `u_wspeed`, `v_wspeed`. Which is exactly what the regrid path
converts -- wind stress, wind velocity and planetary albedo.

So the published path supplies ExoPlaSim's DYNAMICS to EMBM and lets EMBM
compute the heat fluxes with its own tunable transport. The partition stays
where cGENIE can tune it, and the over-transport problem of 28b never arises.

**That gives two architectures, and this project has not chosen between them.**

| | what ExoPlaSim supplies | who computes surface heat flux | partition |
| --- | --- | --- | --- |
| published path | wind stress, winds, albedo, topography | EMBM | tunable in EMBM |
| OCN-10 as written | the full flux set, including net and penetrative shortwave, longwave, sensible and latent heat | ExoPlaSim | not tunable anywhere |

OCN-10 specifies the second: a contract carrying "net and penetrative shortwave,
longwave, sensible and latent heat". That is not what the published path does,
and it is the one that needs OCN-5's loop.

The first has its own cost, which is why this is a fork rather than an obvious
choice. Running EMBM means running a second atmosphere, with `diffamp`,
`diffwid`, `difflin`, `betaz` and `betam` calibrated for Earth, on a planet with
30-hour rotation and 32 degree obliquity. Those would have to be re-tuned -- and
tuning an atmospheric transport to make a coupled answer come out is exactly
what `docs/src/practice/failure-modes.md` class 16 forbids.

So: architecture one avoids the loop and buys a tuning problem this project's
conventions will not permit it to solve the usual way. Architecture two keeps
ExoPlaSim as the only atmosphere and needs the loop. The second is more
consistent with everything else here, but the choice should be made explicitly
rather than by OCN-10 quietly specifying it.

### 28d. And the multiplier is `scf`

`initialise_embm.F:1195` uses `scf` in the wind-stress non-dimensionalisation as
`... *rh0sc*dsc*usc*fsc/(rhoair*cd*scf)`, beside the drag coefficient. It is in
BOTH `ini_embm_nml` and `ini_gold_nml`, and `initialise_goldstein.F:2110` passes
it on as `go_scf`. The shipped configurations set `ea_11` and `go_13` to the
same value, 1.531013488769531300.

So OCN-17's multiplier is the WIND STRESS SCALING, it must be identical in both
components because both consume the stress, and the regrid README's 2.0 and 2.6
are that factor set for an ExoPlaSim-derived stress product against a default
tuned for cGENIE's own. That supersedes the C-grid staggering hypothesis of
section 10c and the gas-transfer coupling of 18c as the primary explanation:
both of those are consequences of the same quantity rather than separate
mechanisms.

### 28e. So the loop is not stylistic, for architecture two

The only resolution is to let the atmosphere respond. Run the ocean, take its
heat transport, return it to the NEXT ExoPlaSim run as a heat-flux convergence
field, and iterate until the partition is self-consistent. That return path
already exists and is what OCN-2 verifies: `nfluko = 1` makes `oceanmod.f90`
read surface code 903 as a monthly W/m2 field and `addfc` apply it under the
header comment "prescribed advection".

OCN-5 declares an iterative loop with an exit predicate rather than a one-shot
forcing, and asserts that this is necessary. This is why: a single pass cannot
converge the atmosphere-ocean transport partition, because the atmosphere in the
forcing climatology has not been told the ocean exists.

It also explains why OCN-2 is ordered before OCN-3 in the row's own reasoning.
The channel that carries the correction back is the thing the loop is built
around, and verifying it with a known field is cheaper than discovering it is
wrong after a circulation exists to blame.


## 29. Ocean albedo: three branches, and the two zenith forms are a free bracket

*Read 2026-08-22, starting from cGENIE's `ocean_alb.F` and following it into
`radmod.f90`.*

### 29a. cGENIE integrates Briegleb over the daylight period

`genie-embm/src/fortran/ocean_alb.F` computes a flux-weighted daily mean ocean
albedo by adaptive quadrature over the daylight period, storing `albo(j,istep)`
per latitude row and timestep. Its integrand `rad_out` is Briegleb:

    rspec = 0.026/(czsol**1.7 + 0.065) + 0.15*(czsol-0.1)*(czsol-0.5)*(czsol-1.0)
    rtot  = rspec + 0.06

So ocean albedo there is a function of solar zenith angle, integrated over the
day. Not of water properties.

### 29b. ExoPlaSim already does this, and ships the same formula as an option

`radmod.f90:2895-2904` sets the direct-beam surface albedo as a sum over surface
types, and the open-ocean term has THREE branches selected by `necham` and
`necham6`:

| switch | open-ocean albedo |
| --- | --- |
| `necham=1`, the DEFAULT | `min(0.05/(zmu0+0.15), 0.15)`, ECHAM-3 |
| `necham6=1` | `0.026/(zmu0**1.7+0.065) + 0.15*(zmu0-1)*(zmu0-0.5)*(zmu0-0.1) + 0.0082`, Briegleb |
| both 0 | the constant `dsalb`, i.e. `doceanalb` |

**This overrides `seamod`.** `seamod.f90:155-158` does apply two namelist
scalars everywhere `dls < 0.5` -- which is how OCN-7 states it -- but `radmod`
then overwrites the open-ocean part. In the default configuration those scalars
survive only for the land and sea-ice fractions, and `doceanalb` contributes
nothing to open water at all: it is a declaration the shortwave replaces. The
sea-ice pair `dicealbmx` and `dicealbmn` does survive, so PHYS-14's sea-ice half
is unaffected. Section 11c carries this where a reader of the flat-array list
will meet it.

Worth noting that the zenith branches are spectrally flat too -- the expressions
for `dsalb(1,:)` and `dsalb(2,:)` are identical -- but here the source SAYS so,
in a comment reading "Currently: we use the same albedo for both spectral
ranges". That is the honest version of failure mode class 31: the comment states
what the code does rather than denying it.

### 29c. The two parameterisations are a free bracket, and it opens where this planet lives

Flux-weighted daily means, computed over the daylight period at both:

| latitude and season | ECHAM-3 (default) | Briegleb | ratio |
| --- | ---: | ---: | ---: |
| equator, equinox | 0.058 | 0.052 | 0.90x |
| 45 deg, equinox | 0.075 | 0.079 | 1.05x |
| 45 deg, solstice at 32 deg obliquity | 0.143 | 0.234 | **1.64x** |
| 60 deg, winter | 0.150 | 0.273 | **1.82x** |
| 70 deg, summer | 0.074 | 0.076 | 1.02x |

They agree within about ten percent wherever the sun is high and diverge to
1.8x where it is low. The cause is the cap: ECHAM-3's `min(..., 0.15)` binds
above 79.5 degrees zenith and Briegleb has no ceiling.

**That divergence region is larger on this planet than on Earth.** At 32 degrees
obliquity against Earth's 23.4, high latitudes spend more of the year in the
low-sun regime where the two disagree, and the winter hemisphere is where sea
ice forms.

So OCN-7's spatial half has a zero-cost A/B available: two parameterisations
already in the source, both defensible, selected by a namelist switch, differing
by up to 1.8x exactly where this planet's obliquity puts more of its year. That
is a better instrument than the row's own framing suggests, and it is a bracket
rather than a fix -- neither branch is obviously right for a K dwarf, and the
question of which is a separate one from whether the term varies at all.


## 30. `genie-plasim`: the coupling contract, written down, and how it is afforded

*Read 2026-08-22. cGENIE ships a PlaSim. It is the same model this project
runs, so the difference between the two trees is the coupling layer and nothing
else -- which makes it a direct readout of what coupling an ocean to this
atmosphere actually requires.*

### 30a. It is UPSTREAM PlaSim, and the coupling surface is 71 lines

`vendor/cgenie/genie-plasim/src/fortran/` is a strict SUBSET of
`vendor/exoplasim/exoplasim/plasim/src/`, plus exactly one file. Absent from it:
`p_exo.f90`, `p_mars.f90`, `p_earth.f90`, the aerosol core, `hurricanemod.f90`,
`glaciermod.f90`, `newsnow.f90`, `carbonmod.f90`, the alternative rain schemes,
and every MPI module. So ExoPlaSim's exoplanet capability IS those additions,
and cGENIE's copy is the Earth model they were added to.

The one addition on cGENIE's side is `geniemod.f90`, 71 lines, and it is pure
declaration -- no logic. **The entire atmosphere-ocean coupling surface is one
array-declaration module.** That is a far smaller thing than the OCN rows have
been assuming.

`lsgmod.f90` on our side is a 4-line stub, not an ocean. Nothing uses it and
nothing can; the name is the only content.

### 30b. The contract, enumerated

Because the module is declaration-only, the contract can simply be read off.

**Ocean to atmosphere, 7 declared and 6 used:** `genie_sst`, `genie_icet`,
`genie_hght_sic`, `genie_frac_sic`, `genie_alb_sic`, `genie_co2`, and
`genie_dflux` marked "not used".

Note the fifth: **sea-ice albedo is supplied BY the ocean model**, so in this
arrangement PlaSim's own `dicealbmx`/`dicealbmn` do not set it. That is the pair
section 29 found surviving `radmod`'s open-ocean overwrite, and PHYS-14 owns it.

**Atmosphere to ocean, 17 fields**, each stored as a full seasonal cycle:
latent and sensible transfer COEFFICIENTS, `netsolar`, `solfor`, `insolar`,
`inlong`, `sat`, `spec_hum`, `pressure`, `evap`, `precip`, `runoff`, wind stress
as `stressx2/y2` and `stressx3/y3`, and `windspeed`. Plus `sfxatm_lnd` for ENTS
carbon.

This is a fuller handoff than the offline regrid path OCN-17 tracks, which moves
wind stress, winds and albedo. The two published paths are different
arrangements, not one; whether EMBM still runs under this one is not settled by
the declaration module and I did not confirm it.

### 30c. The gearing, which is how a coupled ocean is afforded

`plasim.f90:604-655`. With `ngear = 1`, PlaSim integrates for
`ngear_years_plasim` years while accumulating daily means into the seasonal
arrays. For the next `ngear_multiple - 1` blocks it **returns immediately** and
does not integrate at all; the ocean is driven from the stored cycle. The
atmosphere runs one block in `ngear_multiple`.

The non-obvious part, and the reason this is worth recording rather than
inventing later: **the terms that depend on ocean temperature are recomputed
live** against `tstar_ocn` at every geared step -- saturation specific humidity
and hence latent heat, net longwave as `0.98 * sigma * T^4`, and sensible from
`sat_plas - tstar_ocn`. Only the transfer coefficients and the downward
radiative and moisture fields are replayed. The source comment says why:
"Needed for stability." Naive replay of net heat flux would sever the negative
feedback that holds SST, and the scheme is built specifically to keep it.

Evaporation is deliberately NOT recomputed, commented as "necessary for moisture
conservation", with the alternative of rescaling precipitation and runoff
considered and rejected as no better.

### 30d. What a Vesper port would hit

- `NLAT_ATM = 32` is a `parameter`, so the module is compile-time fixed to T21.
  This project runs T42. Mechanical, but it is a recompile and not a namelist.
- `tstar_ocn` is in Celsius here -- `(tstar_ocn + 273.15)**4` -- against
  ExoPlaSim's Kelvin. Rule 3's class of defect, on units rather than longitude.
- Surface emissivity 0.98 and sigma are hardcoded in the geared longwave, so the
  geared branch and the live branch could disagree if either is configured.
- The seasonal arrays are dimensioned `(:,:,360)` by day-of-year. **This is NOT
  a blocker here**: at a 182.8-day orbit and 30-hour rotation this world has
  about 145 solar days per orbit, so the array is oversized and safe. It would
  have to move for a longer year, and it is the same calendar-port class the
  LPJ-GUESS work already went through.


## 31. VOLC-9's weathering columns, confirmed from use -- and the bracket is a different one

*Read 2026-08-22. VOLC-9 required its column reading be confirmed before any
bracket was built on it, on the grounds that a bracket over a misread column is
worse than none. The Fortran that consumes the table settles it.*

### 31a. The columns, from use rather than position

`rokgem_data.f90:814` reads `weath_consts(i,j), j=1,7` positionally, with no
header. What each column does, from where it appears in `rokgem_box.f90`:

| col | meaning | evidence |
| --- | --- | --- |
| 1 | base weathering rate | multiplied into `conv_factor(k)` at 1559, 1592, 1716 |
| 2 | runoff EXPONENT | `loc_runoff(i,j) ** weath_consts(k,2)` at 1568, 1601 |
| 3 | fCa, carbonate fraction | partitions `dum_calcium_flux` at 1848, 1900 |
| 4 | fSi, silicate fraction | partitions `dum_calcium_flux` at 1858, 1910 |
| 5 | osmium yield | 1927-1928 |
| 6, 7 | 187Os and 188Os | 1930-1931 |

The source names columns 3 and 4 outright at 1842: "extra array terms for fCa
and fSi in weath_consts array". **VOLC-9's inferred reading of those two columns
was correct.** The precondition is met.

### 31b. But the conclusion that depended on it does not fire

VOLC-9 expected the two schemes to disagree on fCa/fSi, and reasoned that if so
the disagreement would be the thermostat itself. They do not disagree. The two
files the code actually reads carry IDENTICAL columns 3 and 4:

    carbonate 0.93/0.07   shale 0.39/0.61   sand 0.48/0.52   basalt 0.00/1.00

The 1.00/0.00 values the row cites are in `Amiotte_2003_consts.dat` and
`Gibbs_1999_consts.dat`, which have four columns and which the seven-column
reader cannot consume.

Because there are only TWO selectable 2D schemes, and the source names them:
`rokgem.f90:69` labels `GKWM` as "Gibbs et al (1999)" and `:76` labels
`GEM_CO2` as "Amiotte-Suchet et al (2003)". The four-column files are the
archival originals of the same two schemes -- `GKWM_consts.dat` shares columns 1
and 2 with `Gibbs_1999_consts.dat` exactly, and `GEM_CO2_consts.dat` with
`Amiotte_2003_consts.dat` exactly, apart from one edited acid-volcanics rate,
0.222 against 0.272. So the row's "three published alternatives" are two schemes
in two vintages each.

Note also that the base rates are not comparable across schemes -- 32000 against
1.586 -- because `conv_GKWM` and `conv_GEM_CO2` carry different unit systems.

### 31c. The bracket that is really there: the runoff exponent

The schemes differ in the FUNCTIONAL FORM of the runoff response, which is a
structural difference rather than a constant one:

- **GKWM**, `rokgem_box.f90:1568`: flux proportional to `runoff ** k`, with the
  per-lithology exponents 0.91 carbonate, 0.68 shale, 0.74 sandstone, 0.69
  basalt, 0.75 granite. **Sublinear.**
- **GEM_CO2**, `:1723`: flux proportional to `loc_runoff(i,j)`. **Strictly
  linear**, and its column 2 is uniformly 1.0, consistent and unused.

Both normalise through `r_avg_runoff`, so they agree at average runoff by
construction and diverge away from it. Weighting the exponents by fSi gives an
effective silicate exponent of **0.720**.

Spatially, GKWM against GEM_CO2 silicate flux:

| runoff / mean | 0.2 | 0.5 | 1.0 | 2.0 | 3.0 | 5.0 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| ratio | 1.57x | 1.21x | 1.00x | 0.82x | 0.74x | 0.64x |

And as a thermostat, under a warming-driven runoff increase of factor f:

| f | 1.05 | 1.10 | 1.25 | 1.50 |
| --- | ---: | ---: | ---: | ---: |
| GKWM weaker by | 29% | 29% | 30% | 32% |

**About 30 percent, and flat across the range**, which makes it usable as a
single bracket number rather than a curve. Silicate weathering is the long-term
sink, so this is the thermostat strength directly -- which is what VOLC-9 asked
for, arrived at through the runoff exponent instead of through fCa/fSi.

It bites harder here than on Earth for a specific reason: the divergence lives
entirely in the TAILS of the runoff distribution, since the two forms are pinned
together at the mean. A world whose runoff distribution is shaped differently
from Earth's therefore sees a different scheme-choice penalty, and the shape is
a derived property of the accepted carve and the resulting hydrography rather
than anything readable off the present export.

### 31d. One trap, currently inert

`sub_load_weath` reads seven values list-directed. Given a four-column file it
would continue onto the following line and desynchronise the whole table
silently rather than failing. The archival files are not selectable through
`par_weathopt`, so nothing can trigger this today -- but VOLC-9 currently treats
those files as the scheme definitions, and anyone acting on that would hit it.


## 32. Two one-bucket land surfaces that disagree by an order of magnitude

*Read 2026-08-22. LSHY-3 is specifying a replacement for ExoPlaSim's scalar
bucket and LSHY-4 the vegetation-water loop around it. cGENIE's ENTS is a peer
implementation of the same abstraction, so the two can be compared directly.*

### 32a. The same equation, a different limiter

Both models compute land evaporation as a bulk aerodynamic flux scaled by a
soil-moisture factor. They differ in that factor.

**ExoPlaSim**, `landmod.f90:418`, with `drhsfull = 0.4` at `:53`:

    drhs = min(1, W / (drhsfull * Wmax))

Linear, and **saturated at the potential rate for any bucket above 40 percent
full**.

**ENTS**, `surflux.F:1334`:

    beta = min(1, (W / bcap)**4)

Quartic, reaching the potential rate only at a full bucket.

At equal fractional fill:

| W/Wmax | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.8 | 1.0 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| ExoPlaSim | 0.25 | 0.50 | 0.75 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| ENTS | 0.0001 | 0.0016 | 0.0081 | 0.026 | 0.063 | 0.130 | 0.410 | 1.00 |
| ratio | 2500x | 312x | 93x | 39x | 16x | 8x | 2x | 1x |

They agree only at a full bucket and diverge by an order of magnitude or more
across the dry half of the range. **ExoPlaSim is the permissive extreme.**

Two honest limits on that table. It compares the SHAPE at equal fractional fill,
not absolute evaporation, because the capacity differs too -- ExoPlaSim's `wsmax`
is a uniform namelist scalar while ENTS's `bcap = min(k8, k9 + k10*Csoil)` is
carbon-dependent. And ENTS's `k8`, `k9`, `k10` arrive through the `ents_control`
namelist and are NOT shipped in this tree, so the capacity half of the comparison
cannot be closed from here. Neither limiter is right; this is a bracket, and the
point is that the project's own side sits at one end of it.

The direction matters for more than evaporation: at a given bucket state
ExoPlaSim converts more of its precipitation to evaporation and less to runoff,
and `scripts/error_budget.py` already records that the E/R ratio is amplified
into the carve.

### 32b. Vegetation reaches the ENTS climate through three prescribed-here fields

ENTS carries two prognostic carbon pools, `Cveg` and `Csoil`, and derives three
surface properties from them that ExoPlaSim prescribes:

| property | ENTS | ExoPlaSim |
| --- | --- | --- |
| bucket capacity | `min(k8, k9 + k10*Csoil)`, `setup_ents.F:474` | `wsmax`, uniform namelist scalar |
| roughness | `max(0.001, kz0*Cveg)`, `:476` | `dz0clim`, prescribed field |
| soil albedo | sand-to-peat interpolation in `Csoil`, `initialise_ents.F:239` | `dalbclim`, prescribed field |

Worth being precise about what this does and does not show. It is NOT that ENTS
has a stomatal scheme and ExoPlaSim lacks one -- ENTS has no stomatal term in its
evaporation either. `fv` and `fws` exist but feed photosynthesis, and the
comment at `surflux.F:1386` distinguishing bare-soil evaporation from
evapotranspiration sits above a formula that is the same in both cases. The
whole vegetation-to-water coupling in ENTS runs through `bcap`.

That is useful precisely because it is minimal. `notes/audits/missed-couplings.md`
finding 4 and the error budget both record the one-way coupling as structural,
and LSHY-4 is scoped to close the loop; ENTS shows that an EMIC-class closure of
it is three algebraic functions of two carbon pools, not a land-surface model.


## 33. Egea's family is the shape vocabulary -- but it limits a different flux

*Read 2026-08-22, following section 32, and checked against the paper itself
rather than against the port that cites it -- which is what settles the scope
question below.*

`src/standalone/Vegetation/soil_moisture_stress.jl` offers three models:
`NoMoistureStressModel` with beta = 1; `TuzetMoistureStressModel`, a sigmoid in
leaf water potential needing plant hydraulics; and
`PiecewiseMoistureStressModel`, cited to **Egea et al. (2011)**:

    beta = min(1, max((theta - theta_low)/(theta_high - theta_low), 0) ** c)

### 33a. It is not the same beta

That form is algebraically identical to section 32's two limiters, and it is
tempting to call ExoPlaSim `c = 1` and ENTS `c = 4` in one family. **They limit
different fluxes**, and the paper says so in its own title: "water stress in
coupled photosynthesis-stomatal conductance models".

- ClimaLand's `betam` multiplies PHOTOSYNTHESIS. `photosynthesis_farquhar.jl:335`
  and `pmodel.jl:450` apply it to the carbon rate, and
  `stomatalconductance.jl:141` notes it is "applied to `An` already, so it is not
  applied again here". It is evaluated on ROOT-ZONE AVERAGED soil moisture,
  weighted by root distribution over `canopy.biomass.rooting_depth`,
  `soil_canopy_root_interactions.jl:185-193`.
- ExoPlaSim's `drhs` and ENTS's `beta` multiply the SURFACE EVAPORATION FLUX
  directly, over vegetated and bare ground alike, from a surface store with no
  root weighting.

So the three parameters are a good shape vocabulary for a limiter and the axes
are the right axes, but Egea et al. is not a citation for the bare-soil
evaporation case and must not be used as one.

### 33b. The structural point the comparison actually supports

ClimaLand needs TWO limiters where ExoPlaSim and ENTS have one. Soil water
supply is handled in the soil hydraulics; the canopy carries its own moisture
stress on root-zone moisture; and the two are separate objects. The single
lumped surface beta in the other two models is doing both jobs at once, on a
store with no depth.

That is exactly what `scripts/error_budget.py` records as structural -- "no
stomatal or LAI control and no rooting depth" -- and it now has a concrete
reference implementation to be priced against rather than a description.

### 33c. The third axis, which survives the correction

`theta_low` is named the wilting point or residual water fraction, and ClimaLand
requires `theta_high > theta_low` to sit above the soil's residual water
content. ExoPlaSim and ENTS both have it at zero in their own limiters: each
approaches zero flux only as its store approaches empty, so water held below the
wilting point is still available to be removed. That holds whichever flux the
limiter is attached to, and it is a real axis for LSHY-3 independent of the
scope confusion above.

## 34. Anoxia: computed from air-filled porosity, or switched on by latitude

*Read 2026-08-22. ClimaLand's soil biogeochemistry and the LPJ-GUESS fork solve
the same problem -- when does waterlogging stop decomposition -- by different
kinds of thing.*

### 34a. ClimaLand computes it

`Soil/Biogeochemistry/co2_parameterizations.jl` implements DAMM, the Dual
Arrhenius Michaelis-Menten model of Davidson et al. (2012):

    R = Vmax(T) * MM_sx * MM_o2
    Sx       = p_sx * Csom * D_liq * theta_l**3          ! soluble substrate
    MM_sx    = Sx / (kM_sx + Sx)
    O2_avail = D_oa * O2_f * theta_a**(4/3)              ! Millington-Quirk tortuosity
    MM_o2    = O2_avail / (kM_o2 + O2_avail)

Decomposition is limited at BOTH ends of the moisture range by two separate
mechanisms: substrate diffusion as `theta_l**3` when dry, oxygen supply as
`theta_a**(4/3)` when wet. As air-filled porosity goes to zero, `MM_o2` goes to
zero and respiration stops. **Anoxia is a computed consequence of water
content.**

Worth one note on style: the Arrhenius term is written in centered form,
`Vmax = V_ref * exp(-Ea/R * (1/T - 1/T_ref))`, and the docstring says why --
so that the rate and the temperature sensitivity are approximately orthogonal
under calibration. Algebraically identical, better conditioned to fit.

### 34b. The LPJ-GUESS fork switches it

`somdynam.cpp:502` carries a humped empirical curve in water-filled pore space,
Friend et al. (1997) Eqn 53 after Parton et al. (1993):

| WFPS % | 10 | 30 | 50 | 60 | 70 | 80 | 90 | 100 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| modifier | 0.044 | 0.325 | 0.882 | 0.978 | 0.712 | 0.520 | 0.403 | 0.360 |

It peaks at 60 percent and **bottoms out at 0.360 at full saturation**. A
completely waterlogged mineral soil still decomposes at over a third of its
optimum rate. On this curve alone, saturation does not preserve carbon.

What preserves carbon is a categorical override, and the category is chosen by
latitude. `framework/guess.cpp:1110` and its counterpart partition the single
`PEATLAND` landcover class at `PEATLAND_WETLAND_LATITUDE_LIMIT = 40.0`:

- `is_highlatitude_peatland_stand()`, `lat >= 40`: `moist_mod = RMOIST = 0.4`,
  blending toward `RMOIST_ANAEROBIC = 0.025` once soil carbon passes an acrotelm
  limit of 7.5 kgC/m2 after Wania et al. (2009b). It ALSO sets
  `moist_mod_saturated`, which `somdynam.cpp:599` applies to the slow and
  passive pools so their effective modifier is exactly 0.025.
- `is_true_wetland_stand()`, `lat < 40`: `moist_mod = moisture_modifier(100) =
  0.36` for every pool. `moist_mod_saturated` stays at its default 1.0, so the
  slow and passive pools **never receive the anaerobic treatment**.

### 34c. What that is worth, and why it is a Vesper problem

On the long-lived carbon pools the two branches differ by **0.36 against 0.025,
a factor of 14.4**, and the selector between them is a hard 40-degree latitude
line inside the PEATLAND class. The fast pools are close at low carbon, 0.4
against 0.36, and diverge as peat accumulates -- so the step falls precisely on
the pools that decide whether peat accumulates at all.

A latitude constant is an Earth climatology proxy standing in for where cold,
wet, slowly-decomposing ground occurs. This world has 32 degrees of obliquity, a
182.8-day orbit and a K-dwarf spectrum, so the latitude-to-climate mapping that
makes 40 degrees meaningful on Earth does not carry over. BIO-29 already
requires an inventory of every direct latitude branch and WET-2 already forbids
latitude selection in the wetness classification; this is one such branch with
its consequence priced.

The DAMM contrast says what the replacement looks like: an O2 limitation
continuous in air-filled porosity needs no category and no latitude, and it is
two Michaelis-Menten factors over quantities a soil column already carries.


## 35. Two accelerators in one codebase, and only one of them is guarded

*Read 2026-08-22. **Section 26 already covers what GEMlite IS** -- offline-tracer
acceleration that leaves the circulation alone, its cycling idiom, and the
validity limit its own header states. This section is a different axis on the
same device: section 30c found PlaSim-GOLDSTEIN gearing, and the difference
between how the two accelerators are CONTROLLED is the reusable part.*

### 35a. One thing section 26 did not record

Beyond the validity limit already noted there, `gemlite.f90`'s header warns that
its `ocn` array has module scope only -- "there is an entirely seperate `ocn`
for BIOGEM". Two live copies of the same ocean state, which is the
duplicate-state hazard rule 5 exists for, and it would need a provenance answer
before adoption.

### 35b. The guard is two-sided, and the thresholds are named

`genie-main/genie.F`, with defaults from
`src/xml-config/xml/definition.xml:241-282`:

- **Entry, on a RATE.** `:621-638`. During the full-model phase, if the
  year-on-year pCO2 change exceeds `gem_adapt_dpCO2dt`, default **0.1 ppm/yr**,
  the phase counter is DECREMENTED and the full model runs another year. The
  accelerator does not engage while the system is moving quickly. The config
  description is explicit: "rate-of-change threshold for switching into GEMlite".
- **Extension, on ACCUMULATED ERROR.** `:808-823`. Inside the accelerated phase,
  if cumulative drift from the phase's starting pCO2 stays under
  `gem_adapt_DpCO2`, default **0.1 ppm**, the phase may be extended.
- **A hard cap by default.** That extension happens only if
  `gem_adapt_auto_unlimitedGEM` is set, which is `.false.`, and the source says
  why in place: "some risk of spurious stuff happening if staying in the GEM
  phase for ever". `gem_yr_min = 10` and `gem_yr_max = 1000` bound the phase.

So: engage on a rate, extend on an error integral, and refuse to run unbounded
unless told to, with the reason recorded beside the switch.

### 35c. The contrast, which is the point

The PlaSim-GOLDSTEIN gearing of section 30c has **no adaptive guard at all**.
`ngear_multiple` is a fixed integer and the atmosphere is skipped on a fixed
schedule regardless of what the ocean is doing. The same project built both.

This is a second axis over section 26a's four accelerators. That table sorted
them by IDIOM -- how the jump is taken. This sorts the same devices by whether
anything MEASURES the jump's error while it is being taken, and only GEMlite
does.

The difference is not carelessness: GEMlite has a **scalar that summarises the
accelerated subsystem's drift** and gearing does not. pCO2 is one number, it is
already computed, and its rate and integral both mean something. That is what
makes a guard writable.

The transferable rule for this project's own accelerators is therefore a
question to answer before building one: *what single number measures the drift
this acceleration introduces, and is it already computed?* If there is one, the
guard is two thresholds and a counter. If there is not, the accelerator is a
fixed schedule and its cost is a declared structural term rather than a bounded
one. `config/pipeline.yaml` already carries exit predicates per loop, so the
place to put the answer exists.


## 36. A second TOPMODEL, and where its parameters come from

*Read 2026-08-22. GW-26 takes the saturated-fraction form from PALADYN.
ClimaLand implements it independently, which both confirms the form and shows
the parameter provenance GW-26 has to decide.*

`src/standalone/Soil/Runoff/Runoff.jl` implements SIMTOP, Niu et al. (2005), "A
simple TOPMODEL-based runoff parameterization (SIMTOP) for use in global climate
models":

    f_sat = min(f_max * exp(-f_over/2 * z_wt), 1)          ! :374
    R_ss  = R_sb * exp(-f_over * z_wt)                     ! :428

### 36a. Four things GW-26 does not currently have

1. **A factor of two, and one decay constant for both fluxes.** GW-26 states the
   form as `f_sat = f_sat_max * exp(-f_grad * z_grad)`. SIMTOP's surface decay is
   HALF its subsurface decay -- the same `f_over` appears in both, with `/2` in
   the surface term only. So a single calibrated constant sets both fluxes and
   they are not independent. Bracketing it moves surface and subsurface runoff
   together.
2. **The cap.** `min(..., 1)` is explicit; a saturated fraction can otherwise
   exceed one.
3. **Parameter provenance, split cleanly.** `f_over` and `R_sb` are both
   labelled "calibrated" in their docstrings. `f_max` is not: it is "computed
   from the topographic index CDF per grid cell". That is exactly the split
   GW-26 wants -- the terrain-derived part is computable here from the 15.19 km
   mesh, and the Earth-calibrated part is two scalars, which supports the row's
   existing instruction to DECLARE and bracket `f_grad` rather than fit it.
4. **Ice.** `:246-262` computes the saturated column twice: once including ice
   (`theta_l + theta_i`) for the SURFACE saturated fraction, once with liquid
   only for the SUBSURFACE flux. Frozen ground generates saturation-excess
   surface runoff but not subsurface flow. Neither PALADYN's form as GW-26
   states it nor the row itself carries that distinction, and it matters on a
   world with a cold winter hemisphere.

### 36b. The estimator decision GW-26 flagged now has two peers on the same side

GW-26 records an open decision: PALADYN estimates grid-cell mean water table
from column water content rather than from a lateral solve, and whether that is
reconciled with this project's Dupuit-Forchheimer solution is unsettled.
ClimaLand does the same thing -- `z_wt` is `depth - h` where `h` is the column
integral of a saturation indicator weighted by `(theta - theta_r)/(nu -
theta_r)`, `:246-256`. Two independent implementations, same choice. That does
not settle the decision, but it means the column estimator is the field's normal
practice and the lateral solve is the departure needing the argument.

## 37. Snow: a third density model, and an albedo predictor we do not have

### 37a. Density -- three rungs, and the middle one is cheap

GRAV-8 records ExoPlaSim's constant 330 kg/m3 against PALADYN's prognostic
Kojima (1967) self-loading compaction, which is linear in `g`. ClimaLand sits
between them. `MinimumDensityModel`, `snow_parameterizations.jl:621`:

    rho_snow = rho_min * (1 - q_l) + rho_liq * q_l
    z_snow   = rho_liq * S / rho_snow

Density interpolates from a dry-snow minimum toward liquid water by LIQUID MASS
FRACTION. No load, no compaction, and therefore **no gravity term** -- so it does
not address GRAV-8's argument, which stands. What it does show is that a density
which responds to melt state is available without solving the compaction
problem, and ExoPlaSim's constant misses that effect as well as the gravity one.
ClimaLand's `AbstractDensityModel` is abstract and its docstring anticipates
prognostic variants; only the minimum-density one ships.

### 37b. Albedo -- the two models pick disjoint predictors

**ExoPlaSim**, `landmod.f90:394-405`, is a linear ramp in SURFACE TEMPERATURE
between 263.16 K and `tmelt`, blended by forest fraction:

    zdalb    = (zalbmax - zalbmin) * (dts - 263.16)/(tmelt - 263.16)
    zalbsnow = max(zalbmin, min(zalbmax, zalbmax - zdalb))

**ClimaLand**, `snow_parameterizations.jl:49-52`:

    alpha = min(1 - beta*(rho_snow/rho_liq - x0), 1) * (alpha_0 + d_alpha*exp(-k*cos_z))

Zenith angle and snow density, the latter named "a proxy for grain size and
liquid water content". No temperature term. The two models proxy the same
physics -- grain metamorphism and wet snow -- through different observables, and
one of them is one ExoPlaSim does not have for this surface at all.

### 37c. What the missing zenith term is worth

With ClimaLand's calibrated values, `toml/default_parameters.toml:155-183`:
`alpha_0 = 0.59`, `d_alpha = 0.40`, `k = 1.96`, `beta = 0.97`, `x0 = 0.2`. At
ExoPlaSim's constant 330 kg/m3 the density factor is 0.874.

| cos(zenith) | 1.00 | 0.50 | 0.20 | 0.05 | 0.00 |
| --- | ---: | ---: | ---: | ---: | ---: |
| albedo | 0.565 | 0.647 | 0.752 | 0.833 | 0.865 |

**The zenith term spans 0.30 at fixed density. ExoPlaSim's ENTIRE temperature
range, `dsnowalbmn = 0.4` to `dsnowalbmx = 0.8`, spans 0.40.** The predictor
this model omits carries roughly as much albedo variation as the one it uses.

Flux-weighted daily means at 32 degrees obliquity:

| | equator equinox | 45 equinox | 60 summer | 60 winter | 75 summer |
| --- | ---: | ---: | ---: | ---: | ---: |
| albedo | 0.599 | 0.640 | 0.624 | 0.783 | 0.638 |

A 0.18 seasonal-latitudinal swing that is absent here, and it is largest in the
winter hemisphere at high latitude -- which is where snow is.

Two connections. It lands in the SAME low-sun regime where section 29 found the
two ocean-albedo branches diverging by 1.8x, so ExoPlaSim applies a zenith
correction to open ocean by default while giving snow and ice none: **the same
physical effect is present for one surface and absent for another inside one
model.** And the density factor multiplies the whole albedo, so a wrong density
constant reaches albedo DIRECTLY in peer practice, not only through cover depth
as GRAV-8 currently has it.

The values are Earth calibrations under a solar spectrum and do not transfer to
a K dwarf. The FORM and the magnitude argument do, and both are orthogonal to
the spectral reweighting PHYS-14 already did -- this is a separate axis, not a
correction to that one.


## 38. The leaf-to-canopy bound is asserted, and it is not one-signed

*Read 2026-08-22. ClimaLand's PFT table carries a quantity this project's leaf
dataset does not, and it turns out to decide the sign of an assumption
`analysis/vegetation_albedo.py` currently states without argument.*

### 38a. What the PFT table is

`src/standalone/Vegetation/pfts.jl` is a covarying trait registry of the kind
PCAR-5 and WET-5 describe, and worth noting for its shape alone: a
`pft_param_list` naming every parameter a PFT must define to be valid, each
value carrying its literature source inline, with the sources listed at the top
of the file. Sixteen parameters spanning canopy radiative transfer, conductance,
photosynthesis and plant hydraulics, `rooting_depth` among them.

Its leaf optics are FOUR numbers per class, not two: reflectance AND
TRANSMITTANCE in each of PAR and NIR. CLM5.0 Table 2.3.1, three distinct sets:

| class | a_PAR | a_NIR | t_PAR | t_NIR | t_NIR/a_NIR |
| --- | ---: | ---: | ---: | ---: | ---: |
| needleleaf | 0.07 | 0.35 | 0.05 | 0.10 | 0.286 |
| broadleaf tree | 0.10 | 0.45 | 0.05 | 0.25 | 0.556 |
| grass / crop | 0.11 | 0.35 | 0.05 | 0.34 | **0.971** |

Grass transmits nearly as much near-infrared as it reflects.

As a cross-check, this project's leaf band pair from 553 ECOSTRESS spectra,
`k25v_band1 = 0.1244` and `k25v_band2 = 0.3738`, sits comfortably inside the
reflectance columns. The levels agree; the transmittances have no counterpart
here.

### 38b. Why the missing column decides a sign

`analysis/vegetation_albedo.py:33-35` states: "A canopy treatment would give a
smaller ratio than a leaf one, so the leaf ratio UPPER-bounds the correction",
and the bracket runs from no correction to the leaf ratio on that basis. No
argument is given for the direction.

Canopy albedo depends on the SINGLE-SCATTERING ALBEDO `omega = alpha + tau`, not
on reflectance alone. A leaf that transmits strongly in the near-infrared raises
`omega_NIR` without raising `omega_PAR`, so the canopy retains more of its NIR
albedo relative to PAR than the leaf does -- and since a K dwarf moves flux INTO
the near-infrared, that AMPLIFIES the correction at canopy level instead of
diluting it.

Taking this project's own band pair and CLM's transmittance ratios, through the
semi-infinite isotropic-scattering similarity result
`alpha_canopy = (1 - sqrt(1-omega))/(1 + sqrt(1-omega))`:

| class | canopy k25v/sun | against leaf 1.1501 |
| --- | ---: | --- |
| needleleaf | 1.1379 | smaller, the stated bound holds |
| broadleaf tree | 1.1848 | **larger, the bound fails** |
| grass / crop | 1.2226 | **larger, the bound fails** |

**The assertion holds for one of CLM's three classes and fails for the other
two**, including grass -- which is the class where this project's own per-class
table already shows the largest K-star shift, 0.303 to 0.333.

### 38c. What this does and does not establish

It does NOT give a canopy ratio. The similarity result is a sensitivity
argument, not a canopy model: no leaf angle distribution, no clumping despite
`Omega` sitting right there in the table, no soil background, no finite LAI and
no direct/diffuse split. ClimaLand's actual `TwoStreamModel` has all of those,
and the transmittances are an Earth calibration besides.

What it establishes is narrower and enough: **the sign of the leaf-to-canopy
correction is not determined by leaf reflectance alone**, so an unqualified
"a canopy treatment would give a smaller ratio" is not supported, and the upper
end of the bracket may not be an upper end. Settling it needs leaf
transmittance, which a reflectance library does not carry -- a data gap, not a
derivation gap.

This lands on BIO-18, which is deriving tree and grass endmembers separately.
The class-dependence above is precisely along that split, which makes the
separate treatment more necessary than the row currently argues.


## 39. Soil albedo does not know whether the soil is wet, and the blocker is the state variable

*Read 2026-08-22, closing the ClimaLand sweep.*

`src/standalone/Soil/soil_albedo.jl` offers two parameterizations. The second,
`CLMTwoBandSoilAlbedo`, after Lawrence and Chase (2007) and modified per
Braghiere et al. (2023):

    alpha_band = alpha_band_dry * (1 - S_e) + alpha_band_wet * S_e

with `S_e` the effective saturation averaged over `albedo_calc_top_thickness`,
default **0.02 m**.

Three things follow.

**The band structure already matches.** ExoPlaSim carries `dalbclim1` for below
0.75 um and `dalbclim2` for above, `landmod.f90:118-119` -- the same split this
project uses everywhere. CLM's parameterization is defined on exactly that pair,
so it would drop into the existing structure rather than requiring a new one.

**What is missing is the moisture dependence.** `dalbclim` and its two band
companions are prescribed static fields. Nothing in this project's albedo
products carries a wet/dry axis either -- `analysis/rock_albedo.py`,
`analysis/playa_albedo.py` and the surface-class builder are all moisture-blind.
Wet ground is darker than dry ground, and on land with a strong seasonal wetting
cycle that is a real seasonal albedo term with no representation at all.

**The blocker is the state variable, not the parameterization.** CLM reads
effective saturation over the top two centimetres. Section 32 established that
ExoPlaSim's land is one scalar bucket with no depth, so there is no surface
layer to read and no way to evaluate `S_e` even if the coefficients were in
hand. **A moisture-dependent soil albedo is unavailable here because the state
does not exist, not because the physics is unknown.**

That connects to LSHY-3 beyond its own scope: the replacement land column would
supply exactly the near-surface saturation this needs, so it enables a climate
feedback rather than only improving a hydrological one.

**And it is the second consumer of a state DUST-17 already specifies.** Dust
emission is NOT missing its moisture term -- `aeolian/scripts/build_dust.py`
gates the threshold friction velocity through Fecan et al. (1999) equations 14
and 15, and `aeolian/config/dust.yaml:108-120` already carries the depth and
bulk-density pair that converts a water depth to gravimetric percent, stated
once and reaching the model as the namelist constant `dustwcv`. What DUST-17
owns is the same problem this section describes: the scalar bucket has no
profile, so the emitting layer's water is inferred rather than held. Fecan's
correction and CLM's albedo saturation are both surface-skin quantities and need
not read at the same depth, but they must come off ONE profile -- DUST-17's own
requirement is to avoid "a second central hydrology", and a saturation defined
independently for albedo would be exactly that.

Soil formation is not in this set. `pedology/scripts/build_soil.py:301-303` takes
RUNOFF rather than soil moisture for weathering and argues the choice: leaching
requires water to drain through the profile, and rain that falls and evaporates
carries nothing. That is a deliberate position, not an omission.

No magnitude is quoted because this tree ships none -- the dry and wet values
are spatially varying MODIS-consistent fields rather than constants. Lawrence
and Chase (2007) and Braghiere et al. (2023) are where they come from.


## 40. OCN-12's inventory, geochemistry side: BIOGEM and SEDGEM

*Read 2026-08-22, closing the cGENIE sweep. OCN-12 already has the circulation
side, where GOLDSTEIN's non-dimensionalisation hardcodes Earth's radius and
gravity. The geochemistry side splits differently, and the split is the useful
part.*

### 40a. The calendar is compile-time, and it is stated more than once

`genie-main/src/fortran/cmngem/gem_cmn.f90:585-587`:

    REAL,PARAMETER::conv_yr_d  = 365.25 !360.00 !365.0
    REAL,PARAMETER::conv_yr_hr = 24.0 * conv_yr_d
    REAL,PARAMETER::conv_yr_s  = 3600.0 * conv_yr_hr

Three things in three lines.

**All `PARAMETER`, so compile-time.** Every per-year to per-second conversion in
BIOGEM and SEDGEM runs through `conv_yr_s`, and none of it is reachable from a
namelist. Changing the year length is a recompile, not a configuration.

**`conv_yr_d` carries two commented alternatives on its own line**, 360.00 and
365.0. That is not idle: section 30 found the PlaSim coupling's seasonal arrays
dimensioned `(:,:,360)`, a 360-day year. So the codebase holds a 365.25-day
geochemistry calendar beside a 360-day atmosphere coupling, and the alternative
sits commented out where someone switched it.

**The 24-hour day is hardcoded** as the literal `24.0`, not derived from a
rotation period. This world's is 30 hours.

For scale: `conv_yr_s` is 31,557,600 s against this world's 15,794,006 s, so
every rate constant carried through it is **2.00x** out.

### 40b. A second, inconsistent statement of the same quantity

`sedgem_box_archer1991_sedflx.f90` hardcodes the literal `3.15e7` **twelve
times** as its own seconds-per-year, converting the organic matter decay
constant `rc` -- `par_sed_archer1991_rc`, order 2e-9 per second -- to a per-year
reaction rate.

`3.15e7` is 31,500,000 s. The model's own `conv_yr_s` is 31,557,600. **The
sediment diagenesis module runs on a year 0.183 percent shorter than the rest of
the model.** On Earth that is negligible and invisible. It matters here because
it means the year length is not stated once, so a port that fixes `conv_yr_d`
leaves this module on an Earth year and nothing fails.

`gem_cmn.f90:802` likewise carries `const_rEarth = 6.37E+06` as a `PARAMETER`,
with `6.371E+06` commented beside it -- the same value GOLDSTEIN's `rsc`
hardcodes independently, so Earth's radius is stated at least twice across the
model.

### 40c. The good news, and it changes where the audit's effort goes

BIOGEM's stoichiometry is **not** hardcoded. `biogem_lib.f90:179-184` declares
`par_bio_red_POP_PON`, `par_bio_red_POP_POC`, `par_bio_red_POP_PO2` and
`par_bio_red_PON_ALK` and puts all four in `ini_biogem_nml`, with
`par_bio_red_PC_flex` available to switch C:P to flexible stoichiometry
entirely. OCN-12's "Redfield or fixed stoichiometry" item is already exposed
here, and the fixed ratios are a default rather than an assumption.

**So the porting difficulty is not uniform across this model.** The
biogeochemical parameters are namelist-exposed and the physical scaling
constants are compile-time `PARAMETER`s. An audit that samples uniformly will
spend its effort in the wrong place: the biology is configurable and the physics
core is where the Earth is welded in.


### 40d. "Rate time base" splits three ways, and only one is a unit conversion

OCN-12 asks for an inventory of rate time bases. BIOGEM's fall into three
classes that need different treatment, and the distinction is the deliverable --
a day inside a UNIT is invisible to any search for constants.

**One: physical and kinetic, unambiguous.** `par_bio_remin_sinkingrate` and its
three siblings in `m d-1`, `par_bio_remin_CH4rate` in `d-1`, `par_scav_Fe_Ks`,
`par_bio_remin_opal_K`. These are settling speeds and reaction rates. The "day"
is 86400 seconds and nothing else, so **a port to a 30-hour rotation must NOT
rescale them.** That is a decisive answer for most of the list.

**Two: biological, ambiguous, and a modelling decision rather than a
conversion.** `par_bio_tau`, "biological production time-scale (days)
(OCMIP-2)". Its value was calibrated where the light cycle is 24 hours and
phytoplankton physiology is entrained to it; the model itself resolves no
diurnal cycle, so mechanically the day is again 86400 s. Whether the number
should move with a 30-hour rotation cannot be settled by unit analysis and
should be recorded as an open choice rather than silently converted.

**Three: not temporal at all -- spectral, and one-signed.** Two parameters set
the light response and neither has a time base:

- `par_bio_c0_I`, half-saturation for light in **W m-2**, Doney et al. (2006).
- `par_bio_I_eL`, **light e-folding depth in metres**, OCMIP-2.

Both are solar-spectrum calibrations. `par_bio_c0_I` is the marine twin of the
constant BIO-25 owns on land, where LPJ-GUESS's `CQ = 4.6e-6` is documented "for
solar radiation at 550 nm": a light response fixed on an irradiance whose photon
content and PAR fraction both change under a K dwarf.

`par_bio_I_eL` is the sharper of the two, because its error has a known sign.
It is a SINGLE SCALAR attenuation depth, and water's absorption is strongly
wavelength dependent -- transparent in the blue-green, strongly absorbing in the
red and near-infrared. A K dwarf puts more of its flux where water absorbs more,
so **the photic zone is shallower here than one solar-calibrated e-folding depth
implies, and marine production is confined nearer the surface.** The direction
follows from the spectrum alone. **The magnitude is computable from what is
already in the tree, and it is large.**

### 40e. What the redder star costs the photic zone: 0.71

Section 12b already argued this direction qualitatively -- "water absorbs
strongly beyond about 600 nm, so the light that reaches a phytoplankton cell is
bluer than the light at the surface". What follows is its magnitude, and the
distinction between this figure and section 12a's 0.869 is drawn there.

`exoplasim/data/water/hale_querry_1973_liquid_water.dat` is this project's own
extraction of Hale and Querry (1973) Table I, deliberately covering 0.75 to
4.0 um. Converting to an absorption coefficient, `a = 4*pi*k/lambda`:

| lambda um | 0.75 | 0.80 | 0.90 | 1.00 | 1.20 | 1.50 | 2.00 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| e-folding depth m | 0.383 | 0.509 | 0.147 | 0.028 | 0.010 | 0.0008 | 0.0001 |

**The longest e-folding depth anywhere above 0.75 um is 0.51 m.** Every photon
in shortwave band 2 is deposited in the top metre or two of water as heat. What
reaches the photic zone is band 1 and essentially nothing else, so the band-1
flux fraction IS the fraction available below the surface layer:

| | band-1 fraction |
| --- | ---: |
| k25v, real spectrum via `lib/stellar.band_fractions` | **0.382** |
| Sun, 5772 K Planck -- `analysis/vegetation_albedo.py`'s own reference | **0.537** |
| ratio | **0.712** |

The k25v figure reproduces the project's canonical `flux_fraction_band1 = 0.382`
exactly, so this uses no new number.

**This star delivers 71 percent as much light to the water below the first metre
or two as the Sun does**, before any question of how the remaining visible
attenuates with depth.

Three limits on that figure, none of which change its sign. It is the flux
reaching below the surface layer, NOT the photic depth itself -- the depth also
needs visible attenuation, and the visible half of Table I is not extracted here
(the existing file stops at 0.75 um deliberately, and the PDF is a scan whose
table would not extract reliably, so it was not guessed at). It is pure water,
whereas attenuation in productive water is dominated by chlorophyll and CDOM,
which absorb more strongly in the blue and therefore partly offset the penalty.
And it treats band 2 as wholly removed, which the table above justifies.

**It is also not only a biology result.** The same calculation says a larger
share of this star's shortwave is deposited as heat in the topmost water. For
the present slab ocean that is invisible, since the slab is one layer. For
OCN-1's multi-layer column it is not: penetrative solar heating would be more
surface-concentrated here than Earth tuning assumes, which stratifies the top of
the column more strongly. OCN-10 already names "net and penetrative" solar in
its forcing contract, so the place to carry it exists.
