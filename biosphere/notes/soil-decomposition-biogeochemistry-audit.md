# Soil decomposition and biogeochemistry audit

**Scope:** the active CENTURY C-N-P path in the vendored LPJ-GUESS-CNP fork,
from litter entry through soil-organic-matter turnover, mineral nutrients,
plant uptake and the pedology feedback. This was a static source and literature
audit. No build, LPJ-GUESS execution or other material computation was run.

Fire transfer is covered by `fire-model-audit.md`; plant hydraulics, regolith
depth and groundwater exchange are covered by
`plant-hydraulics-groundwater-audit.md`. They appear here only where they set a
boundary condition for decomposition or nutrient cycling.

## Active model path

The Vesper instruction path enables `ifcentury 1`, nitrogen limitation and soil
N transformations, while phosphorus limitation remains off pending BIO-5
through BIO-10. The CENTURY path in `modules/somdynam.cpp` carries eleven bulk
pools: surface and soil structural and metabolic litter, fine and coarse woody
debris, surface and soil microbial pools, surface humus, and slow and passive
SOM. First-order decay is modified by temperature, upper-soil moisture,
lignin, texture and nutrient availability. Flexible pool C:N and C:P ratios
couple decomposition to immobilisation and mineralisation.

This is a defensible baseline model class, not a generic description of soil.
Parton et al. (1993) developed and tested the source architecture against
Earth grasslands and 0--20 cm soil observations. The named microbial pools are
conceptual turnover pools: their size does not control enzyme activity or
substrate consumption as explicit microbial biomass does.

After a sampled forcing window, `equilsom()` replaces daily integration with a
40,000-year monthly accelerator. That accelerator is part of the active C-N
baseline, not only the future P experiment. Its equations must therefore be
equivalent to the daily operator before any equilibrated result is accepted.

## Findings

### 1. The equilibrium accelerator does not preserve the daily nutrient operator

The daily path stores mean daily uptake and leaching fractions. The monthly
accelerator then applies each mean once rather than composing the daily
survival fractions. More seriously, mineral N and P leaching subtract
`mass * (1 - mean_leach_fraction)`: a small daily leaching fraction therefore
removes nearly the entire monthly pool instead of a small compounded fraction.
The same aggregation mismatch affects uptake.

Phosphorus forcing is sampled incorrectly as well. `soilpadd()` adds the
year-to-date accumulators `apwtr` and `apdep` to their purported means every
day. This forms a triangular sum of cumulative totals rather than an annual
input. `equilsom()` divides weathering by the number of sampled years, but does
not divide `apdep_mean` at all. Inside the accelerated solve, P uptake,
leaching and additions mutate `pmass_labile` directly, bypassing the
`pmass_add()` labile--sorbed equilibrium used by the daily path.

These are implementation-correctness defects, not reasons to prefer a more
complex soil model. SDEC-1 requires a deterministic daily-to-monthly operator
fixture, C-N-P closure and parity before patching individual expressions. It is
a blocker for the present C-N result as well as C-N-P activation.

### 2. The phosphorus topology was an unbounded drain, and is now the published one

The daily path computed a reversible-looking transfer
`USORB * sorbed - USSORB * strongly_sorbed`, removed that amount through
`pmass_add()`, and reported it as an external `P_SOIL` flux. It never
incremented `pmass_strongly_sorbed`, which is declared, initialised, serialised
and reported but was assigned nowhere in the tree. With that stock pinned at
zero the reverse term vanished, so the line was not an exchange between two
pools: it was a first-order drain of the sorbed pool that could never shut off.

The magnitude, from the `Spmax` of 77 to 145 gP/m2 the fork adopts from Wang et
al. (2010) Table A1: 0.52 to 0.97 gP/m2 per Earth year, against a default
`soiltype.pwtr` weathering supply of 0.003 gP/m2 per Earth year. Under
`ifplim 1` the soil phosphorus system could not have reached a steady state.

It did not surface as a conservation failure because `Patch::pcont()` excluded
the destination pool, so the phosphorus had genuinely left the accounted system
and the books balanced, and because `MassBalance::check_patch_P` is declared in
`guess.h` and called from nowhere in this tree.

The topology chosen is the published one: Dantas de Paula et al. (2025) Appendix
A2 and Wang et al. (2010) Eq. D10, a reversible exchange between the sorbed and
strongly sorbed stocks. `pmass_strongly_sorbed` is assigned, so both terms are
real and the stock fills to the size of the sorbed pool, where the net flux goes
to zero. The transfer is internal, so it is no longer reported as a soil P loss,
and `Patch::pcont()` counts the pool instead; double-booking it would break the
balance the two together close.

Occlusion is a DECLARED ABSENCE rather than the alternative topology. `UOCC` had
no citation anywhere in the tree, occluded phosphorus is terminal so an
unbracketed rate would set this world's long-run soil phosphorus stock, and
`equilsom` spins the pools for 40000 model years without saving and restoring
the sorbed stocks the way it does `pmass_labile`, so a terminal sink inside that
device would drain every cell over a span the run does not represent. `UOCC` is
deleted and `Patch::pcont()` counts `pmass_occluded` anyway, so that adding a
flux later is a change to one file. The argument is in
`biosphere/notes/phosphorus-cycle-parameterisation.md`.

`USORB` and `USSORB` are 0.0067 per EARTH year, and now divide by
`VESPER_EARTH_YEAR_DAYS` rather than by the model year, under
`biosphere/notes/time-base-unit-contract.md`. Their equality is Wang et al.
(2010)'s own and is not a copying artifact. `PMASS_SAT` is Parton, Stewart and
Cole (1988) Fig. 3 exactly and is not a copy, but the labile P pool it gates is
the fork's wider Hedley-labile one and not that figure's resin-extractable
orthophosphate, so its ramp saturates everywhere. `PCONC_SAT` has no phosphorus
source at all and its ramp can never be reached. Both are their own rows.

### 3. A bulk SOM column breaks the depth, root and groundwater contracts

LPJ-GUESS has fifteen standard soil layers, but CENTURY SOM is not indexed by
depth. All root litter enters the same `SOILSTRUCT` and `SOILMETA` pools. One
temperature near 25 cm and upper-soil water status control the whole mineral
soil column. Natural mineral soil does not receive a water-table or oxygen
state, so decomposition and N transformations cannot respond to saturated
deep roots, capillary access, anoxia or redox boundaries except through the
separate wetland/peatland model forms.

This conflicts directly with Vesper's 0.02--5 m regolith field and with
PLHY-3/PLHY-4's pending root-depth and groundwater exchange. Koven et al.
(2013) show why layered hydrology plus bulk biogeochemistry is not a neutral
abstraction: vertical litter input, transport, oxygen and depth-dependent
protection change both stocks and their response. Their Earth implementation
is evidence for the missing coordinate, not a drop-in parameter set.

The pedology feedback currently makes the asymmetry worse. It reads the bulk
`cpool.out` soil C total, then LPJ's organic-soil-property path redistributes
that total over its physical layers with the fixed Earth-global
`beta_SOC = 0.976` profile and dumps the residual into the deepest layer. That
profile changes hydraulic and thermal properties but does not become the
biogeochemical profile that generated the carbon. SDEC-3 therefore owns a
vertically resolved model-form bracket tied to the same regolith/root/water
geometry as PLHY-3 through PLHY-5; it must not carry the Jobbagy--Jackson curve
forward as an observation of Vesper.

### 4. Mineral protection is reduced to Earth texture regressions

Transfer from microbial and slow pools to slow/passive SOM is controlled by
fixed linear functions of clay or clay plus silt inherited from CENTURY. The
pedology model already distinguishes andic material and phosphate fixation,
but `SoilInput` skips the `andic` and `pfixation` columns. It also carries no
iron/aluminium oxide, allophane, aggregate-capacity or polyvalent-cation term
into carbon stabilization.

Cotrufo et al. (2013) and Lehmann and Kleber (2015) both place microbial
products, mineral association, aggregation and accessibility at the centre of
stable SOM formation. Mineralogy can make soils with similar texture behave
differently; Vesper's explicitly generated volcanic/andic contrast is exactly
such a case. SDEC-4 links the pedology mineral state to both P sorption and an
explicit mineral-protection bracket. This extends BIO-5/BIO-6 rather than
creating a second soil-map interface.

### 5. Litter chemistry is fixed globally and only lignin:N controls its split

`somdynam.h` fixes lignin fractions at 0.20 for leaves, 0.16 for roots and 0.30
for wood for every PFT; the wood value still carries a source TODO. The
structural/metabolic split uses lignin:N, while the source itself asks whether
an equivalent P control is needed. Root chemistry has no depth or absorptive
versus structural-root distinction.

These are Earth plant proxies layered beneath the already-declared use of Earth
PFTs. They matter more once GRAV-7 changes wood and structural-root investment,
and they also control fire fuel quality. SDEC-5 makes tissue-specific litter
chemistry part of the PFT trait registry and brackets the fixed-CENTURY and
trait-resolved cases without inventing a uniquely Vesperian chemistry.

### 6. Nutrient limitation lacks active acquisition processes

In the fork, low mineral N or P can reduce donor-pool decomposition when
immobilisation demand cannot be met. Plants otherwise take up the labile
mineral pools. There is no explicit mycorrhizal investment, phosphatase or
organic-acid release, root exudation/priming trade-off, or biochemical
mineralisation of organic P. Dantas de Paula et al. identify these omissions as
a likely source of excessive progressive P limitation; Wang et al. include an
explicit biochemical P-mineralisation route based on acquisition cost.

SDEC-6 registers a no-retuning acquisition sensitivity coupled to root carbon
cost, mineralogy and microbial demand. It is a model-form bracket after the
central P pool is correct, not a hidden tuning knob for BIO-9's prediction.

### 7. Microbial physiology can reverse the soil-carbon response

The CENTURY microbial pools participate in fixed transfer fractions and
stoichiometric constraints, but microbial biomass does not catalyse decay and
there is no carbon-use-efficiency, enzyme or microbial-turnover response.
Wieder et al. (2013) demonstrate that an explicit microbial operator can change
the sign or magnitude of soil-carbon responses to warming and extra litter;
they also show that growth-efficiency temperature response is itself poorly
constrained. Their implementation omits C-N interactions, mineral protection
and coarse woody debris, so it cannot simply replace CENTURY here.

SDEC-7 therefore keeps CENTURY as the auditable baseline and defines one
explicit-microbial/mineral-stabilisation structural sensitivity, with microbial
growth efficiency and adaptation registered as priors. It is uncertainty on
soil/climate feedback, not a prerequisite for proving the adapter runs.

### 8. Soil N transformations assume an Earth gas and redox environment

The active `ntransform.cpp` path partitions nitrification and denitrification
with fixed Earth-calibrated temperature, water-filled-pore-space and pH
response curves and constants from `global_soiln.ins`. It uses upper-soil water
status but no atmospheric pressure/O2 boundary, gas diffusivity, water-table
redox state or depth structure. Pedology pH reaches the code, but Vesper's
stronger-gravity pressure field currently does not.

Because gaseous loss changes mineral N availability, this is not merely an
emissions-reporting issue. SDEC-8 joins BIO-23's pressure/O2 handoff and
SDEC-3/PLHY's saturation state to a registered N-transformation parameter and
model-form bracket. Methane remains the separate, disabled BIO-28 scope.

### 9. Current artifacts cannot diagnose or close this system

`biosphere/scripts/run_lpj_guess.py` retains total `cpool.out`, `nsources.out`
and `cflux.out`, but not the per-pool C-N-P stocks, litter chemistry,
heterotrophic respiration, mineralisation/immobilisation, P sorption and sink,
leaching, soil state, N gases or accelerator-versus-daily residual needed to
locate a budget error. Warning-only internal balance checks are not an
acceptance artifact.

SDEC-9 extends BIO-10/BIO-14's acceptance contract and the BIO-15
soil--biosphere convergence assessor. It must distinguish real convergence of
the coupled pedology--biosphere iteration from convergence onto the imposed
Earth vertical profile or a fast-solver artifact.

## Order of work

1. Fix and fixture the accelerator under SDEC-1 before accepting even the C-N
   baseline.
2. Reconcile P topology under SDEC-2, then complete BIO-5 through BIO-10.
3. Add SDEC-9 diagnostics alongside those correctness fixes so closure is an
   artifact rather than a console warning.
4. Close SDEC-3/SDEC-4/SDEC-8 interfaces with the existing pedology,
   groundwater and pressure work before interpreting deep-root, groundwater-fed
   or nutrient-loss behavior.
5. Treat SDEC-5 through SDEC-7 as pre-registered model-form sensitivities after
   a correct baseline exists.

No task above authorizes an LPJ-GUESS run. Source work and no-simulation
fixtures can proceed while the machine is occupied; any model execution still
requires explicit permission.

## Primary sources read

- Parton et al. (1993), *Observations and Modeling of Biomass and Soil Organic
  Matter Dynamics for the Grassland Biome Worldwide*,
  `10.1029/93GB02042`.
- Wang, Law and Pak (2010), *A global model of carbon, nitrogen and phosphorus
  cycles for the terrestrial biosphere*, `10.5194/bg-7-2261-2010`.
- Cotrufo et al. (2013), *The Microbial Efficiency-Matrix Stabilization (MEMS)
  framework integrates plant litter decomposition with soil organic matter
  stabilization: do labile plant inputs form stable soil organic matter?*,
  `10.1111/gcb.12113`.
- Wieder, Bonan and Allison (2013), *Global soil carbon projections are
  improved by modelling microbial processes*, `10.1038/nclimate1951`.
- Koven et al. (2013), *The effect of vertically resolved soil biogeochemistry
  and alternate soil C and N models on C dynamics of CLM4*,
  `10.5194/bg-10-7109-2013`.
- Lehmann and Kleber (2015), *The contentious nature of soil organic matter*,
  `10.1038/nature16069`.
- Dantas de Paula et al. (2025), *Including the phosphorus cycle into the
  LPJ-GUESS dynamic global vegetation model (v4.1, r10994) - global patterns
  and temporal trends of N and P primary production limitation*,
  `10.5194/gmd-18-2249-2025`.
