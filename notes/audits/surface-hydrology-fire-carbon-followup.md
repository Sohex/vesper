# Surface, hydrology, fire and carbon follow-up audit

This note audits five proposed omissions against the repository and the
vendored model sources.  It concerns the fictional Vesper simulation.  No
ExoPlaSim, hydrography, dust or LPJ-GUESS integration was run.

## Sources and acquisition status

The source audit covered ExoPlaSim's snow-albedo calculation, both dust
emission paths, the BLAZE effects code and surface builders, the lake
representation note and carve code, the weathering/outgassing diagnostic, and
Orogen's exported tectonic fields.

The following primary papers were added to `references/`, read in full on
2026-08-21 and indexed in `references/INDEX.md`:

- Betts and Ball (1997), *Albedo over the boreal forest*,
  `10.1029/96JD03876`;
- Essery (2013), *Large-scale simulations of snow albedo masking by forests*,
  `10.1002/grl.51008`;
- Randerson et al. (2006), *The impact of boreal forest fire on climate
  warming*, `10.1126/science.1132075`;
- Jin et al. (2012), *The influence of burn severity on postfire vegetation
  recovery and albedo change during early succession in North American boreal
  forests*, `10.1029/2011JG001886`;
- Bernus and Ottlé (2022), *Modeling subgrid lake energy balance in ORCHIDEE
  terrestrial scheme using the FLake lake model*,
  `10.5194/gmd-15-4275-2022`;
- Marty and Tolstikhin (1998), *CO2 fluxes from mid-ocean ridges, arcs and
  plumes*, `10.1016/S0009-2541(97)00145-9`; and
- Le Voyer et al. (2019), *Carbon Fluxes and Primary Magma CO2 Contents Along
  the Global Mid-Ocean Ridge System*, `10.1029/2018GC007630`.

The papers constrain mechanisms and the shape of appropriate tests.  Their
Earth vegetation, fire histories, lake depths and tectonic fluxes are not
Vesper parameters.

## Findings

### 1. Forest snow masking exists, but its structural and spectral port does not

The proposed missing pathway is already active.  Surface code 212 supplies
`dforest`.  With the production configuration's `NVEG=0`, `landmod.f90` blends
forested and unforested snow-albedo endpoints by that fraction before blending
snow with the background surface.  `build_surface_albedo.py` writes code 212
from modelled tree cover and already forces lakes and barren fractions out of
the canopy.  A uniform snow layer therefore does **not** simply repaint all
trees white.

The residual is still material.  The mask is one scalar tree-cover fraction;
it has no canopy exposure relative to snow depth, plant/stem area, sky-view or
gap fraction, intercepted snow, unloading, branch burial, standing dead
structure or PFT distinction.
The two forested snow endpoints are initialized as fixed fractions of the
K-star snow endpoints (`albsmaxf=0.5*albsmax`, with the minimum set from that),
not derived from K-star radiative transfer through the accepted canopy.  Betts
and Ball show the size and vegetation dependence of the Earth phenomenon.
Essery shows that gap fraction, canopy snow and ground snow are separable model
states, but also found that a simple PFT-weighted scheme, a gap-fraction scheme
and a two-stream scheme gave similar large-scale results when supplied with
realistic land cover and parameters.  The first port therefore needs explicit
exposure and snow states plus a hierarchy of reduced operators; a full
two-stream canopy is not automatically the best central model.  Those sources
justify a structural model and tests, not transplanting their Earth albedos.

This belongs in `BIO`, because LPJ canopy structure supplies the state and
ExoPlaSim remains the radiative consumer.  It must extend BIO-17/BIO-18 rather
than replace the existing mask.

### 2. Fire has no persistent surface-disturbance state

BLAZE retains daily and annual burned area and applies combustion, litter and
mortality effects.  It does not retain char cover, bare exposure, standing
dead structure or time since fire as a surface-boundary state.
`build_surface_albedo.py` and `build_surface_roughness.py` reconstruct their
fields from surviving tree/grass cover, substrate, lakes and relief.  FIRE-8
and FIRE-9 name albedo as a consumer, but neither defines the missing state or
its recovery.  This is a real interface gap.

The proposed sign needs correction.  A severe fire can darken a snow-free
surface immediately, but canopy loss exposes snow and can raise spring albedo
for years to decades.  Jin et al. found an initial summer decrease followed by
recovery, while the spring increase appeared immediately and grew with burn
severity; Randerson et al. found the multidecadal albedo term could dominate
the fire-cycle radiative balance.  The response also changes roughness,
interception and evapotranspiration.  A single `char -> bright bare soil`
timescale or a universally warming sign would therefore encode the wrong
model.  The accepted artifact needs severity, char/bare/standing-dead cover and
age, with optical loss and ecological succession separated and all fractions
closing.

### 3. Both dust paths use ExoPlaSim's scalar bucket

The concern is correct.  The offline dust builder reads climatological `mrso`
and converts the water depth to gravimetric moisture using a declared emitting
depth and bulk density before applying the Fecan threshold.  The resident
interactive patch gathers `dwatc` directly and applies the same kind of scalar
conversion.  Neither path sees groundwater-fed surface wetness, an explicit
top layer, liquid/ice partition or a wet/dry subgrid mosaic.

The proposed ownership is slightly wrong.  Dust should not read an independent
LPJ soil column or map `PLHY-4` directly into ExoPlaSim.  LSHY-3 already owns the
planned climate-side land column, LSHY-4 makes it the central vegetation-water
owner, and LSHY-6 preserves lake/playa/upland/groundwater-fed fractions.  The
durable fix is for offline and interactive dust to consume the same
provenance-bound **top emitting-layer liquid-water state** from that column,
including bulk density, texture, frozen water, capillary/groundwater supply and
emitting-area fraction.  The old `mrso`/`dwatc` result remains the exact scalar
reduction, not a second central hydrology.

### 4. Lakes have the wrong thermal state, not literally zero thermal mass

`exoplasim/notes/lake-representation.md` already identifies the absence of a
lake model, water-column temperature and lake ice.  The current subgrid lakes
are land cells with modified albedo and bucket capacity, so they inherit a
soil thermal column; they do not have zero heat capacity or necessarily track
air temperature synchronously.  The problem is that the thermal material,
depth, mixing, phase and surface-flux state are those of land rather than a
distinct lake.  Majority-lake mask flips can use the slab ocean, and the slab's
mixed-layer array is already spatially indexed, but neither route resolves the
dominant subgrid population without tiles.

CLIM-26's roughly -0.9% land-precipitation moisture response did not bound this
error.  It measured an annual moisture term, whereas lake depth controls the
phase and amplitude of surface temperature, evaporation and ice.  Bernus and
Ottlé found seasonal temperature differences of several kelvin for large/deep
lakes and materially different evaporation when a lake energy balance replaced
bare soil.  On Vesper, the correct first step is an offline/reduced seasonal
bound using each basin's area and hypsometric depth, not an immediate global
lake-model port.  Only a decision-relevant bound should trigger a coupled
experiment, and every such run still requires explicit permission.

FLake in Bernus and Ottlé is a freshwater model, and the paper identifies
salinity as one cause of lake-temperature error because it changes albedo,
evaporation and seasonal thermal behaviour.  That matters especially here:
closed-basin Vesper lakes can become saline or briny.  The reduced bound must
therefore cross depth and ice physics with a freshwater-to-brine property
bracket rather than treating the freshwater result as the central answer.

The carve criterion computes Penman open-water evaporation from climatological
atmospheric forcing rather than reading the land bucket.  Thermal inertia can
therefore move the verdict only through seasonally resolved lake evaporation
and coupled changes to precipitation/runoff; it should be tested specifically
on basins close to their recorded evaporation margin, not asserted from lake
area alone.

### 5. The required global outgassing is already diagnosed; tectonic capacity is not

`weathering_fluxes.py` already reports silicate consumption and the implied
steady-state outgassing relative to Earth.  `config-rationale.md` explicitly
treats fixed 450 ppm as a checked assumption, compares its required flux with a
planet-mass heat-supply scaling, and lists spreading rate, mantle temperature,
gravity and ocean fraction as unresolved.  The volcanic-sulfate builder already
consumes the same implied outgassing ratio.  The proposal is therefore not a
missing carbon-balance diagnostic.

A narrower tectonic plausibility check is worthwhile.  Orogen exports boundary
class, plate Euler poles and per-region surface velocity, so ridge/arc length
and relative motion can be measured on the same 10M support.  Its angular
velocity is explicitly in arbitrary units, however, not rad/Myr.  Dividing the
weathering sink by fault length would yield a useful required mol/km/year
diagnostic, but **cannot falsify 450 ppm by itself**.  Carbon flux depends on
magma production or convergence throughput, mantle/slab carbon, volatile
heterogeneity, degassing efficiency, subduction retention/recycling, and
plume/rift sources.  Observed Earth segment fluxes vary by orders of magnitude.
Le Voyer et al. measured more than three orders of magnitude of segment-scale
flux variation; normalizing their tabulated results by segment length still
spans more than two orders of magnitude because primary-magma carbon and magma
production vary independently.  Marty and Tolstikhin likewise derive ridge,
arc and plume fluxes from melt emplacement, source carbon and degassing, with
most arc carbon recycled from the slab.  These are direct reasons not to turn
fault length into capacity.

The honest task is to partition the required global source across ridge, arc,
rift/plume and metamorphic brackets; report length and relative-throughput
intensities; and state which absolute-rate conversion remains declared because
the terrain generator has no physical tectonic clock.  It complements rather
than replaces the existing global Earth-ratio check.  If every plausible
source/efficiency bracket fails, fixed CO2 is inconsistent; a large per-km
number alone is not that test.
