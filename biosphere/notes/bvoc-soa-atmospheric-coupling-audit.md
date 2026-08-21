# BVOC, secondary organic aerosol and atmospheric coupling audit

Worldbuilding. Vesper is invented; every atmosphere and biosphere quantity below
is a model input or prediction, not an observation.  This audit was written
2026-08-21 from the vendored LPJ-GUESS-CNP source, the Vesper adapter and
ExoPlaSim aerosol patches, and from the primary papers listed at the end.  No
LPJ-GUESS build or simulation and no climate run was performed.

## Executive finding

LPJ-GUESS can provide a useful terrestrial **BVOC source sensitivity**, but
turning `ifbvoc` from 0 to 1 would not produce a defensible SOA forcing.  Four
different models are presently being collapsed into that one switch:

1. plant production and emission of speciated volatile carbon;
2. oxidation by OH, O3 and NO3 in a NOx-dependent atmosphere;
3. production, partitioning, growth, transport and removal of organic aerosol;
4. direct radiation and aerosol--cloud interactions.

Only the first exists, and even it is an Earth-calibrated diagnostic that is not
closed in the plant carbon budget.  The existing ExoPlaSim work supplies a good
carrier for a **prescribed direct-aerosol sensitivity**: the multi-species patch
has a fourth slot, per-species optical ratios, vertical scale heights and column
fields.  It supplies neither atmospheric chemistry nor aerosol microphysics,
and its shortwave aerosol effect is explicitly omitted from cloudy fractions.

The practical route is therefore staged: correct and retain the emission
inventory; pass it through a reduced, offline chemistry--transport bracket;
route the resulting seasonal aerosol optical-depth field into the fourth
radiative slot; and keep cloud/CCN effects as a separate structural sensitivity.
The atmosphere must be iterated because climate controls vegetation and
emissions, while the resulting aerosol and oxidant changes feed back on climate.

## What exists now

### LPJ-GUESS source

- `vendor/lpj-guess/data/ins/global.ins` sets `ifbvoc 0`.  All annual and
  monthly isoprene/monoterpene output declarations are commented out.
- `biosphere/scripts/run_lpj_guess.py` does not retain BVOC files.
- `biosphere/scripts/build_lpj_driver.py` already carries ExoPlaSim's daily
  temperature range, and `modules/vesperinput.cpp` interpolates it into
  `climate.dtr`.  BVOC leaf temperature is the only current consumer.
- `modules/bvoc.cpp` implements Arneth et al. (2007) isoprene and Schurgers et
  al. (2009) monoterpene production/storage.  Internally it reports isoprene and
  nine monoterpene compounds.  `commonoutput.cpp` collapses those to total,
  endocyclic (`MT1`) and other (`MT2`) outputs.
- The Vesper PFT file supplies Earth emission capacities.  Several entries say
  explicitly that values were copied from another PFT, and every PFT uses the
  same compound-by-compound storage-fraction vector.

### Atmospheric carrier

`exoplasim-3.4.2-multi-species-aerosol.patch` raises the radiative capacity to
four species: dust, sea salt, sulfate and one spare.  A prescribed species can
receive a gridded band-1 column optical depth, a scale height, two-band
single-scattering/backscatter/extinction ratios and a thermal-IR absorption
ratio.  That is enough to conduct a direct-radiation experiment after a burden
and optical properties exist.

It is not an aerosol model.  It does not oxidize VOCs, predict gas/particle
partitioning, form or grow particles, transport separate organic size modes, or
activate cloud droplets.  `vendor/exoplasim/docs/aeromod.rst` also states, and
`radmod.f90` implements, that shortwave aerosol scattering and absorption act
only in the cloud-free fraction.  Longwave extinction is additive with cloud
transmission, but that does not repair the missing cloudy-sky shortwave term.

### Atmospheric composition

`config/planet.yaml` prescribes ozone and fixes CH4 at 1.6 ppmv from Rugheimer
et al. (2013)'s 5000 K case.  That paper holds modern-Earth biogenic surface
fluxes fixed.  The configuration already says this is an assumption pending a
biosphere source.  BVOCs consume OH, alter methane lifetime, and can either form
or destroy tropospheric ozone depending on NOx.  A biosphere-derived BVOC source
combined with fixed-Earth oxidants and fixed-Earth methane is not a closed
atmospheric state.

## Findings

### 1. Activation is an output and acceptance problem, not a switch

The dormant module is mechanically reachable, but the run harness would discard
its files and the standard outputs lose most compound identity.  This matters
before chemistry: alpha-pinene and beta-pinene have different highly oxidized
molecule yields, while isoprene can suppress monoterpene nucleation in some
regimes.  A total-terpene carbon flux cannot reconstruct those pathways.

The activation gate needs compound-resolved monthly fluxes, PFT and patch
provenance, storage-pool state, complete-cell checks and a source-carbon ledger.
It must remain closed until findings 2--4 have at least declared brackets.

### 2. Earth PFT averages hide the dominant source uncertainty

`initbvoc()` turns standard emission capacities into electron fractions at
370 ppm CO2, 30 C, 1000 micromol photons m-2 s-1, 12 h daylight, a fixed first
canopy-layer absorption fraction and the model's Earth photon conversion.
Those are standard measurement conditions, not Vesper plant traits.

Arneth et al. (2007) evaluated the isoprene form at five Earth ecosystems and
emphasized canopy composition and disturbance history.  Schurgers et al. (2011)
then showed why a PFT mean is especially weak: BVOC emission is not a trait by
which the PFTs were constructed, species aggregation over Europe tended to
overestimate isoprene and underestimate monoterpenes, and apparent capacities
varied by up to about 50% and 40% respectively as composition changed.  MEGAN2.1
puts global isoprene uncertainty near a factor two, many monoterpenes near a
factor three, and local/seasonal uncertainty higher still.

Vesper can retain transplanted Earth biochemistry as its central case, but the
emission capacities and storage strategies belong in PCAR-5's covarying live
plant trait registry with emitter/non-emitter and species-composition brackets.
They must not be retuned to obtain a desired aerosol forcing.

### 3. Emitted carbon is reported but never removed from the plant

`bvoc()` calls `Individual::report_flux()` for isoprene and nine monoterpenes.
No corresponding carbon is removed from assimilation, NPP, a labile pool or the
monoterpene storage pool's parent plant carbon stock.  `monstor` changes the
timing of diagnostic emitted carbon but is itself created outside the plant
carbon ledger.  The model can therefore emit carbon without lowering biomass or
raising ecosystem carbon loss.

The process formulation allocates a fraction of photosynthetic electron
transport to terpene synthesis, so energy is also being used.  A coupled run
must charge both the emitted/stored carbon and its declared energetic cost,
without subtracting it once as carbon and again as a generic respiration cost.
The existing `monstor` serialization must be covered by restart-continuity
fixtures and the storage stock must close over daily, orbit and
equilibrium windows.  This extends PCAR-10 rather than bypassing it.

### 4. The environmental response inherits the unresolved physiology port

Several constants are explicitly Earth-specific:

- `daytime_temp()` maps daylength to radians with `/24`, while Vesper rotates in
  30 h and a forcing step covers 0.8 rotation;
- the standard conversion assumes a 12 h light period and Earth `CQ`;
- seasonality uses an 11 h photoperiod threshold, a 5 C threshold, GDD and a
  fixed 0.05 d-1 decay;
- monoterpene storage uses 2, 80 and 365 **day** time constants;
- `leafT()` fixes air density at 1.204 kg m-3 and uses PFT-constant aerodynamic
  conductance; it receives net shortwave and actual ET but not ExoPlaSim's
  downwelling/net longwave, pressure, humidity or wind;
- BVOC production uses non-water-stressed electron transport and actual
  water-stressed assimilation only as an on/off `adtmm > 0` gate.

The first four fall under BIO-21/BIO-22, PCAR-1 and BIO-25.  Leaf energy needs
BIO-23's local pressure, radiation, humidity and wind.  Continuous hydraulic
stress should follow PLHY-1 rather than the binary gate.  The simple inverse
`370/co2` response should remain a named Arneth sensitivity, not silently become
a universal response outside its tested CO2 and plant domain.

### 5. Isoprene plus nine monoterpenes is not the whole volatile source

MEGAN2.1's global inventory is roughly 50% isoprene and 15% monoterpenes; it
represents about 150 compounds and identifies important methanol, acetone,
acetaldehyde, sesquiterpene and stress-induced classes.  Some lower-flux
compounds have much higher aerosol-forming capacity than isoprene.  Heat,
drought, ozone, herbivory, mechanical damage and fire can also change the amount
and composition of emissions, while the retained LPJ module represents almost
none of those induced pathways.

Implementing 150 compounds would add false precision.  The minimum honest
product is a named `LPJ isoprene+MT only` lower model and a lumped missing-VOC
source bracket, split at least into non-SOA oxygenates, sesquiterpene/high-yield
organics and induced stress emissions.  FIRE-6 supplies the damage/fire seam;
PLHY supplies drought state.  These terms must remain separate from primary
fire smoke.

### 6. Vesper oxidants cannot be replaced by a fixed Earth SOA yield

Isoprene and terpenes react with OH by day and with O3/NO3 on other pathways.
NOx changes radical recycling, ozone production and product volatility; aerosol
acidity, sulfate, water and ammonia open further partitioning pathways.
Shrivastava et al. (2017) shows that even the direction of several NOx effects
is pathway-dependent.  Vesper's K-star UV spectrum changes photolysis and OH,
while its current ozone and methane are prescribed outputs of a separate
fixed-Earth-flux photochemical experiment.

The first atmospheric layer should therefore be a reduced oxidant bracket,
driven by Vesper UV/ozone, water vapour, temperature and the soil/fire NOx
sources that the project actually adopts.  It must return both precursor loss
and changes to O3, OH/oxidizing capacity, CO and methane lifetime.  A fixed
mass-yield calculation is useful only as a deliberately broad screening bound,
not as the central coupling.

### 7. SOA burden requires chemistry, transport, partitioning and removal

The 31-model AeroCom comparison found natural+anthropogenic SOA sources spanning
13--121 Tg/yr, OA lifetimes spanning 3.8--9.6 days and more than an order of
magnitude diversity in vertical OA.  Models with semivolatile treatment
produced much larger SOA sources than simpler fixed-yield forms.  It also found
that global models generally underestimated observed OA.

The reduced Vesper operator therefore needs, at minimum: precursor-specific
OH/O3/NO3 loss, a volatility or reversible-partitioning representation, HOM/new
particle production as a separate branch, condensation onto existing dust,
salt and sulfate surface area, horizontal and vertical transport, wet/dry
removal, and carbon conservation between gas, aerosol and deposition.  It can
run offline on archived ExoPlaSim meteorology and emit a seasonal prescribed
column/vertical climatology.  Reusing a dust lifetime or transport field would
be wrong: precursor lifetimes are hours and OA removal/partitioning is
composition and cloud dependent.

### 8. SOA is not generically an absorbing aerosol

The existing CLIM-29 text combines smoke and SOA and says both absorb with a
forcing opposite to sulfate and salt.  That is incorrect.  Shrivastava et al.
(2017) reports that most climate models treat OA as nearly white and that
scattering dominates their direct SOA effect.  Some SOA becomes wavelength-
dependent brown carbon, but absorption depends on precursor, NOx, ammonia,
aqueous processing, aging and photobleaching.  Alpha-pinene+OH SOA is cited as
nearly nonabsorbing, unlike several anthropogenic aromatic pathways.  Both the
sign and magnitude of direct SOA forcing remain low-confidence.

Vesper needs a K-star-weighted optical bracket spanning mostly scattering SOA,
aged/hygroscopic particles and a weakly absorbing biogenic brown-carbon end.
Primary fire black/brown carbon remains FIRE-9/CLIM-29's distinct absorbing
source.  The fourth ExoPlaSim species can carry the direct SOA cases after the
cloudy-sky shortwave omission is either repaired or explicitly bounded.

### 9. Particle number and clouds are a separate structural problem

Mass and optical depth do not determine cloud condensation nuclei.  Kirkby et
al. (2016) demonstrated pure-biogenic nucleation from highly oxidized
alpha-pinene products, including an ion-induced path without sulfuric acid.
Gordon et al. (2016) put that mechanism into an aerosol microphysics model and
found a potentially large change in preindustrial aerosol--cloud forcing.
Carslaw et al. (2013) found 45% of the variance in industrial-period indirect
forcing arose from uncertain natural emissions, including BVOCs.

Those Earth forcing numbers are not Vesper priors.  They establish that omitting
the pathway cannot be disguised as a direct-optics approximation.  ExoPlaSim
has no prognostic aerosol activation or droplet-number response.  The project
needs a separately scored structural bracket: no cloud effect; an offline CCN
bound; and, only if decisive, a minimal aerosol-number/cloud-albedo port.  Size,
hygroscopicity, ion production, seed aerosol and cloud supersaturation must be
carried explicitly enough that an optical-depth multiplier cannot stand in for
them.

### 10. This creates a new convergence loop and accepted artifact

The causal path is circular:

    climate -> vegetation/physiology -> BVOC flux and speciation
            -> oxidants/SOA burden and size -> radiation/clouds -> climate

Every crossing must record the source-climate, biosphere, chemistry and optics
hashes and preserve Vesper's absolute-time/orbit convention.  Convergence must
be assessed on vegetation, precursor flux, aerosol burden/optical depth,
oxidant/trace-gas changes and the climate fields that drive them.  A final
experiment should be a preregistered matched set: no BVOC feedback,
emissions-only, direct-SOA scattering/absorption bracket, and any accepted CCN
branch.  No LPJ-GUESS or climate experiment is authorized by this audit.

## Recommended order

1. Land the source ledger, compound outputs and Earth-PFT uncertainty registry.
2. Correct the inherited physiology/time/pressure dependencies before enabling
   the source in an interpreted run.
3. Build a reduced oxidant and SOA chemistry--transport operator that outputs a
   prescribed seasonal climatology and conservation report.
4. Price direct radiation with the existing fourth aerosol slot and corrected
   optics, keeping primary smoke separate.
5. Use the result to decide whether a cloud/CCN port can change a project
   decision; do not build it solely because Earth models contain one.
6. Close the full loop only after each stage has a no-simulation contract and an
   accepted artifact.

## Sources read

- Arneth et al. (2007), *Process-based estimates of terrestrial ecosystem
  isoprene emissions: incorporating the effects of a direct CO2-isoprene
  interaction*, `10.5194/acp-7-31-2007`.
- Schurgers et al. (2009), *Process-based modelling of biogenic monoterpene
  emissions combining production and release from storage*,
  `10.5194/acp-9-3409-2009`.
- Schurgers et al. (2011), *Effect of climate-driven changes in species
  composition on regional emission capacities of biogenic compounds*,
  `10.1029/2011JD016278`.
- Guenther et al. (2012), *The Model of Emissions of Gases and Aerosols from
  Nature version 2.1 (MEGAN2.1): an extended and updated framework for modeling
  biogenic emissions*, `10.5194/gmd-5-1471-2012`.
- Carslaw et al. (2013), *Large contribution of natural aerosols to uncertainty
  in indirect forcing*, `10.1038/nature12674`.
- Tsigaridis et al. (2014), *The AeroCom evaluation and intercomparison of
  organic aerosol in global models*, `10.5194/acp-14-10845-2014`.
- Kirkby et al. (2016), *Ion-induced nucleation of pure biogenic particles*,
  `10.1038/nature17953`.
- Gordon et al. (2016), *Reduced anthropogenic aerosol radiative forcing caused
  by biogenic new particle formation*, `10.1073/pnas.1602360113`.
- Shrivastava et al. (2017), *Recent advances in understanding secondary
  organic aerosol: Implications for global climate forcing*,
  `10.1002/2016RG000540`.

The PDFs are held under `references/` and recorded as read in
`references/INDEX.md`.  Rugheimer et al. (2013), already held and read there,
supplies the separate fixed-biogenic-flux photochemical state used by the
current atmosphere configuration.
