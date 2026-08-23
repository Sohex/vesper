# Abiotic nutrient delivery audit

**Recorded:** 2026-08-21  
**Scope:** abiotic nutrient movement from rock, sediment and atmosphere to the
root zone used by the vendored LPJ-GUESS-CNP fork, including weathering, dust,
lightning, runoff, groundwater and depositional renewal. Fire recycling is
covered by `fire-model-audit.md`; biological fixation and acquisition remain
biotic processes. This was a static source, artifact and literature audit. No
ExoPlaSim, chemistry or LPJ-GUESS simulation was run.

## Result

The repository does not yet have an abiotic nutrient-delivery model. It has
several useful pieces, but no mass-carrying path joins them:

- pedology assigns rock phosphorus content and a relative release score;
- the dust model emits total mineral mass deposition without composition or
  solubility;
- hydrography routes water topology, and groundwater routes only water;
- the ecological driver supplies one spatially uniform declared nitrogen
  deposition scalar; and
- LPJ-GUESS accepts immediately available mineral N and P inputs, then treats
  leaching as an external loss at each independent grid cell.

Most importantly, `phosphorus_budget.py` does not consume the dust-deposition
artifact despite its pipeline dependency and README description. Its aeolian
result is only the ratio of source-rock P concentration to land-mean rock P.
Its weathering and basin-delivery values are also relative ranks, not
kgP/m2/time. They are useful hypotheses about sign and geography, but cannot
be passed to LPJ-GUESS as either weathering or deposition fluxes.

The central architecture should be one phase- and species-resolved nutrient
ledger. Each source must state material, chemical form, mass per receiving
area per absolute time, interval, source and destination. Dissolved species
follow the accepted land-water ledger; particulate material follows erosion,
transport, deposition and dissolution; retained mineral stocks are distinct
from immediately labile nutrients. LPJ-GUESS should receive only the flux or
stock that actually crosses into its declared soil pools, never a total dust
or rock inventory relabelled as bioavailable input.

## Sources read

The source trace covered pedology weathering, solute-routing and phosphorus
scripts; the aeolian dust model and pipeline contracts; the ecological driver;
hydrography and groundwater artifacts; and the CNP fork's weathering,
deposition, snow storage, mineral pools and leaching paths.

The audit read the repository copies of Hartmann and Moosdorf (2011), Hartmann
et al. (2014), Okin et al. (2004), Mahowald et al. (2008) and Hudson-Edwards et
al. (2014). Four additional primary papers were fetched into `references/`,
read, and recorded in `references/INDEX.md`:

- Porder et al. (2007), on uplift, erosion, soil residence time and root-zone
  phosphorus;
- Aciego et al. (2017), on measured dust and bedrock nutrient supply;
- Schumann and Huntrieser (2007), on lightning NOx production; and
- Ardaseva et al. (2017), on composition-dependent lightning chemistry in
  Earth-like atmospheres.

Their Earth fluxes and calibrations are constraints and model-form evidence,
not measurements of Vesper.

## Findings

### 1. No artifact closes source material to a root-zone nutrient input

`pedology/scripts/phosphorus_budget.py` says that it computes weathering,
routing and aeolian return. In fact it multiplies P concentration in ppm by a
relative lithology release factor. It groups that relative quantity by whether
a cell drains to a closed basin or the ocean and concentrates it geometrically
over basin floors. The only surface-water use is to reject currently wet
basins as dust sources. The script never opens `dust_baseline.nc` and carries
no dust mass, source provenance, particle size, P phase or dissolution.

That makes the current basin result a geomorphic hypothesis rather than a
budget. Endorheic retention does not itself put P in a root zone: material may
remain dissolved in lake water, precipitate in salt, bind to sediment, become
biogenic material, or be buried beneath an inaccessible floor. Likewise, a
source lithology richer or poorer in P than the land mean says nothing about
the mass deposited at a destination.

There is no equivalent abiotic N ledger. Lightning is being investigated for
fire ignition, but no flash-to-fixed-N path exists. The LPJ driver instead
declares one global deposition value with no source or atmospheric budget.
ANUT-1 establishes a common ledger before individual sources are implemented.

### 2. The weathering fields do not determine an absolute P supply or initial stock

Pedology's WHAK-style `weathering_intensity` is explicitly normalized and
dimensionless. Hartmann et al. (2014) supply the source lithology, runoff,
temperature and soil-shielding architecture, while Hartmann and Moosdorf
(2011) calibrated absolute catchment P-release equations on Japan. The latter
reach up to 390 kgP/km2/Earth-year in that domain, but remain Earth catchment
relationships whose runoff, exposed lithology and soil shielding must be
reconstructed rather than copied as Vesper constants.

An absolute rock-to-solution flux also needs a material boundary. Parent P
concentration multiplied by a dimensionless intensity cannot say how much
rock reacts per area and time. Erosion, regolith production, exposed mineral
surface, soil depth and residence time determine whether fresh apatite enters
the weathering zone. Porder et al. show that both rapid erosion, which exports
material before it dissolves, and very long residence, which permits depletion
and occlusion, can lower root-zone available P.

The fork sharpens this omission. `Soil` initializes its labile, sorbed,
strongly-sorbed and occluded P stocks to zero. During spin-up, weathering and
deposition construct the mineral inventory from flux alone. A 40,000-year
accelerated SOM solve is not a pedogenic history, and existing accelerator
defects are separately tracked by SDEC-1. ANUT-2 therefore owns an absolute
weathering/renewal artifact, while ANUT-10 owns the distinct initial-stock and
history interface. BIO-5 consumes ANUT-2 rather than trying to convert the
present relative score directly.

### 3. Total dust mass is not bioavailable atmospheric phosphorus

`aeolian/scripts/build_dust.py` retains size internally for settling, but its
soil-facing output is total deposition in g mineral/m2/Earth-year. It does not
retain source tags or size-resolved deposition, and it supplies no elemental
composition, mineral phase, wet/dry partition or dissolution kinetics.

The distinction is load-bearing. Okin et al. demonstrate that dust-P can
maintain terrestrial ecosystems over long timescales, and Aciego et al. found
measured dust supply equal to modern erosional output and 10--20% of estimated
millennial bedrock-P input at their montane sites. Mahowald et al. estimate
mineral aerosol as 82% of global atmospheric total P, yet their modeled dust
phosphate source is only one tenth of dust total P. The repository's direct
Bodele measurements are more restrictive: Hudson-Edwards et al. found only
2--4% water-soluble P in source sediment and dust, with much of the P in
authigenic apatite.

LPJ's `snow_pinput()` cannot represent that material. It stores all `dpdep` in
snow and then `pmass_add()` places it directly into the labile--sorbed mineral
system. Passing bulk dust P there would make every deposited apatite particle
immediately available. ANUT-3 must produce separate total particulate,
dissolved/bioavailable and retained mineral inputs with explicit uncertainty;
BIO-8 may carry only the accepted receiving-pool flux.

### 4. A lightning diagnostic is not yet a nitrogen fixation flux

FIRE-1's planned convection-based lightning potential is the right upstream
starting point, but fire and nitrogen need different reductions. Ignition
needs cloud-to-ground flashes and a successful-ignition efficiency. Abiotic N
fixation depends on total discharge energy and the yield of NOx per energy,
including intracloud flashes.

Schumann and Huntrieser synthesize an Earth best estimate of about 5 TgN/year
with an uncertainty factor near five and show that flash rate, energy per
flash, NO yield and vertical placement all matter. Ardaseva et al. show why the
yield is not portable as one Earth scalar: the cooling and reaction paths
depend on bulk composition and density, and their N2/O2 and N2/CO2 cases
produce different NO chemistry.

Fixed N must then be oxidized, transported and removed as NOy by wet or dry
deposition. Co-location with convective precipitation matters, as do local
pressure, atmospheric O2/CO2/N2, photochemistry and the absolute-time basis.
ANUT-4 derives a composition-aware total-flash-to-NOy source and deposition
bracket from FIRE-1/EFOR rather than turning the fire ignition field into
fertilizer.

### 5. The LPJ nitrogen boundary invents amount, chemistry, geography and time

`build_lpj_driver.py` defaults to 0.5 kgN/ha/Earth-year and correctly labels it
declared rather than measured. The binary carries that one scalar for every
cell and forcing interval. `vesperinput.cpp` divides it equally between NH4 and
NO3 and applies it through the model calendar; BIO-24 already owns the
orbit-versus-absolute-time correction.

The equal split is not chemically neutral. Lightning initially supplies
oxidized N, not ammonium. A natural NHx field requires an identified reduced-N
source and atmospheric path; biological fixation inside LPJ is a separate
biotic source and cannot be used to justify atmospheric NH4. The stock
LPJ-GUESS input machinery already distinguishes wet/dry NHx and NOy, proving
that the model boundary can carry the necessary species even though the
Vesper adapter does not.

ANUT-5 replaces the scalar with a spatial, chronological, species- and
wet/dry-resolved interface, while retaining zero or explicit brackets for
sources the project has not modelled. It must not tune deposition to compensate
for the empirical Cleveland biological-fixation relation.

### 6. Nutrient losses leave LPJ but do not enter rivers or groundwater

LPJ grid cells are hydrologically independent for nutrients. Mineral NO3 and
labile P leaching are proportional to LPJ percolation and reported as external
losses; organic C and N can also leave, while organic P is presently treated
differently. No exported mass is routed to a downstream cell, floodplain,
lake, groundwater store or ocean, and no dissolved nutrient can return with
baseflow or capillary rise.

Upstream routing cannot close that gap. `solute_routing.py` routes selected
major ions with a local runoff proxy but does not carry N or P. Hydrography and
the groundwater solver carry water only. The phosphorus script separately
routes a relative score, so connecting it to LPJ would double-model a source
while still losing LPJ's actual leachate.

LSHY-2/LSHY-3 first define surface runoff, drainage/recharge and groundwater
exchange without double-counting `P-E`. ANUT-6 then follows dissolved species
with that ledger and particulate nutrients with physical erosion/sediment,
including sorption, residence and accessible deposition on river, floodplain,
lake and basin-floor surfaces. PLHY-4 supplies the lower-boundary water
exchange; nutrient return is a separate transported mass, not an inference
from water-table depth.

### 7. Geological renewal and episodic deposition are outside the current model

Static lithology is not a nutrient-supply history. Uplift, erosion, regolith
production and burial renew or isolate parent material; Porder et al. show that
their effect depends on the ratio of soil residence time to dissolution and
occlusion time. Current pedology explicitly has no time coordinate, and the
phosphorus basin geometry does not say how long enrichment accumulated.

Fresh volcanic ash and tephra are also distinct from the crystalline rock of
the same nominal lithology. The repository already records much faster glass
dissolution and large silica release from fresh ash, but models volcanic
sulfate aerosol only, not nutrient-bearing ash deposition. Sea-salt aerosol is
likewise radiative-only; over land it is primarily a salinity/base-cation
input, not a generic fertilizer. Fire ash is nutrient recycling and remains
FIRE-7's scope.

ANUT-7 registers these as explicit model-form/source screens. It should add a
source only if its plausible mass can change the root-zone ledger, and retain
event timing where pulse delivery matters. It must not infer volcanic nutrient
flux from sulfate optical depth or sea-salt nutrient flux from aerosol
extinction.

### 8. C-N-P limitation does not prove general nutrient sufficiency

The CNP fork can limit plants by N and P only. Pedology already represents or
routes material relevant to K, Ca, Mg, S and other solutes, while dust and ash
can import those elements plus Fe and trace nutrients. None enters a plant
availability pool. This does not justify immediately adding every element to
LPJ-GUESS, but it does make “nutrient-limited biosphere” broader than the
implemented claim.

ANUT-8 performs a conservative adequacy screen for K, S, Ca, Mg, Fe and
biologically important trace elements against parent stocks, plausible
delivery, loss and plant demand. The durable result is a declared CNP model
boundary plus fail-closed brackets for any plausible co-limiter. Ocean nutrient
delivery and marine productivity are outside the terrestrial LPJ result and
must be stated as a separate boundary rather than silently counted as export.

### 9. Acceptance requires stock-and-flux closure across component boundaries

No current artifact can answer how much source P or fixed N entered the
atmosphere/soil boundary, changed phase, remained in soil, entered biomass,
left by leaching/erosion, accumulated downstream or reached the ocean. Relative
ranks and concentration ratios cannot close mass, and a locally closed LPJ
budget would still omit exported nutrient destinations.

ANUT-9 extends BIO-10/BIO-14 and SDEC-9 with source-to-destination checks,
species/phase reconciliation and reduced fixtures. Closure should be tested at
each handoff and over the whole domain, with residuals separated from declared
terminal sinks. These fixtures require no climate or LPJ execution. Any
production response experiment remains subject to explicit permission.

## Order of work

1. Define ANUT-1's ledger and ANUT-9's reduced closure fixtures first.
2. Build absolute weathering under ANUT-2, initial stocks under ANUT-10 and
   atmospheric particulate delivery under ANUT-3 before BIO-5/BIO-8 map their
   fluxes into LPJ.
3. Reuse FIRE-1's upstream convection work for ANUT-4, but derive total-flash
   chemistry independently of FIRE-3's ground-flash ignition conversion.
4. Land ANUT-5's N interface with BIO-24's absolute-time contract.
5. Add ANUT-6 transport after LSHY has assigned each water flux one owner.
6. Use ANUT-7/ANUT-8 as pre-registered source and co-limitation screens, not
   as permission to tune additional nutrients to a desired productivity.
