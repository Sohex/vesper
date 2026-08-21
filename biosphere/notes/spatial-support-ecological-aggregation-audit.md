# Spatial support and ecological aggregation audit

**Recorded:** 2026-08-21  
**Scope:** the spatial chain from the fully resolved 10M-region World Orogen
export through ExoPlaSim, hydrography, groundwater, pedology and the vendored
LPJ-GUESS-CNP fork, including feedback and scoring products. This was a static
source, artifact and literature audit. No World Orogen, ExoPlaSim or
LPJ-GUESS simulation was run.

## Result

T42 is an operational convention, not a physical or architectural limit. The
repository already has terrain exports and climate binaries for the ladder
T21/T42/T85/T127/T170, and optimization is making the upper rungs increasingly
practical. Resolution should therefore be selected by convergence of the
quantities the coupled biosphere consumes and returns, with cost recorded as a
constraint rather than embedded as a scientific rule.

The reference beneath that ladder is `precarve-craton-10m`: about ten million
native regions at a 7.60 km mean edge. The generator's designed terrain
information bottoms out near 20 km, so this mesh oversamples the information
floor by about 2.6 times. Its preserved-basin catalogue is also near its
measured crossover: 3,621, 6,345, 9,419 and 9,649 basins at 2.5M, 5M, 10M and
25M regions, with only 2.4% gained after 10M at much higher cost. This makes the
10M export the correct fine-support reference for aggregation tests. It does
not make T170 a 7.6 km climate model.

The present pipeline loses more information through support semantics than
through the choice of T42 alone. It turns partial land into a binary climate
cell, mixes parent materials before nonlinear pedogenesis, assigns one dominant
LPJ soil category and one climate history to a cell, treats demographic patches
as if they could stand in for spatial response units, and scores extensive
outputs against full binary land-cell area. Several surface operators also
aggregate parameters before applying nonlinear flux laws. A finer grid reduces
some of these errors but does not define the missing operators or preserve the
discarded covariances.

The central correction is a spatial-support contract plus a small set of
conservative, field-specific operators. Every artifact must state its geometry,
area measure, fractions, intensive/extensive/categorical semantics, source and
destination support, and whether a nonlinear process was evaluated before or
after aggregation. Where cell means are insufficient, carry a limited mosaic
of co-located terrain, soil, hydrologic and ecological response units with area
weights. LPJ stochastic patches remain samples of demography within one
environment; they are not those response units.

Finally, Loop D cannot be a cosmetic resolution change after Loops A--C have
settled. A finer climate grid changes coastline, orography, precipitation,
runoff and recharge, then soil, vegetation, albedo, roughness and the carve
verdict. The existing restart-converter specification provides the right
continuation mechanism: map the accepted coarse state into the next support,
then settle Loop A and every affected downstream loop on that support. The
converted state is an efficient initial condition, not evidence that the finer
coupled state is already at equilibrium.

## Support ledger

These are different scales and must not be called simply "the resolution":

1. **Native terrain support:** the ~10M-region Orogen mesh, 7.60 km mean edge,
   which fully samples the generator's ~20 km information floor.
2. **Atmospheric numerical support:** the Gaussian T21/T42/T85/T127/T170
   ladder on which ExoPlaSim solves climate.
3. **Land response support:** partial coast, lake, barren/rootable, terrain,
   soil and hydrologic classes within an atmospheric cell.
4. **Ecological demographic support:** independent LPJ grid cells and the
   stochastic patches/stands sampled inside each one.
5. **Comparison and decision support:** the common grid or regions on which
   equilibria and extensive totals are compared.

For orientation, not as a statement of the shortest dynamically resolved
wavelength:

| climate grid | cells | equatorial longitude-cell width | ~10M regions per global cell, mean |
| --- | ---: | ---: | ---: |
| T21, 32x64 | 2,048 | 751 km | 4,883 |
| T42, 64x128 | 8,192 | 375 km | 1,221 |
| T85, 128x256 | 32,768 | 188 km | 305 |
| T127, 192x384 | 73,728 | 125 km | 136 |
| T170, 256x512 | 131,072 | 94 km | 76 |

The widths use Vesper's 7,645.2 km radius. Gaussian meridional spacing is
comparable away from the poles. The region counts are global means; coast and
land-cell counts differ. Even T170 therefore aggregates many native regions,
and the spectral model's effective dynamical resolution is coarser than one
longitude cell.

## Sources read

The source trace covered the build/grid resolver and mesh binning library;
boundary-condition, roughness and albedo builders; hydrography coupling and
surface-water paths; pedology's lithology mixing, subgrid slope and nonlinear
weathering; the LPJ driver, soil selection and scoring code; the pipeline loop
documents; the resolution/terrain audit; the restart-converter specification;
and the ExoPlaSim binary matrix.

Four primary papers were fetched into `references/`, read and recorded in
`references/INDEX.md`:

- Rastetter et al. (1992), which formalizes ecological aggregation error as
  the difference between applying a nonlinear fine-scale relation to means and
  taking the expectation of the fine-scale response, and develops moment,
  partition and calibration corrections;
- Avissar and Pielke (1989), which evaluates homogeneous subgrid surface
  classes separately and area-weights their atmospheric fluxes;
- Wood et al. (1988), which shows that a representative hydrologic support is
  process- and topography-dependent rather than a universal grid size; and
- Essery et al. (2003), which directly compares tiled and aggregate land
  surfaces and shows that coupled atmospheric feedback can amplify differences
  that appear small offline.

Their Earth grids, classes and response magnitudes are model-form evidence, not
measurements or preferred scales for Vesper.

## Findings

### 1. The repository has a five-rung ladder but an operational T42/T85 ceiling

`source/README.md` and `rebuild_binaries.py` define the intended ladder as
T21/T42/T85/T127/T170. `precarve-craton-10m` carries all five Gaussian exports.
In contrast, pipeline products and consumers still name T42 or T85 literally:

- hydrography defaults to coupling matrices for T42 and T85 only;
- surface water, carve verdict and carve-list export default to T42 coupling;
- pipeline writes for boundary conditions, albedo, roughness, soil water and
  aerosols are rooted under `inputs/t42`;
- configuration independently states `resolution`, `latitudes` and
  `longitudes`, so a valid truncation can be paired with stale dimensions; and
- some mesh consumers open `exoplasim-T42` by name even where they mean "the
  export that carries raw mesh data."

The last convention can remain, because raw mesh data are stored once in the
T42 directory, but it needs to be named as a storage carrier rather than
allowed to imply model resolution. SPAT-2 makes resolution a parameter and one
registry the source of grid geometry.

### 2. The 10M mesh is the reference, but the reusable support artifact is missing

`lib/gridding.py` correctly assigns every native region to a Gaussian cell and
can compute a land fraction, land-weighted mean or land fraction of a class.
Several consumers independently rebuild counts, moments or class shares, and
GRID-2 records that hypsometry is computed three times and persisted nowhere.

The shared artifact should be emitted for every requested ladder grid from the
10M mesh. It needs total and land area, native-region count, partial land,
rootable, barren and solved-lake fractions, elevation moments/hypsometry, and
co-located terrain/lithology/hydrologic class shares. Empty-cell fallbacks must
be recorded and tested rather than silently filled from the nearest region.
At 10M, previous comments about 35 empty T85 polar cells on the old 2.5M mesh
are not current measurements and must not be carried forward.

SPAT-3 extends GRID-2 rather than creating a competing hypsometry product.

### 3. Binary coastline conversion changes area and physics

`build_boundary_conditions.py` first computes partial land fraction, then
thresholds it into a binary ExoPlaSim land mask. Elevation is averaged only over
the land portion and assigned to the entire accepted land cell. The operation
therefore changes coastline, effective land area, surface reservoir ownership
and coastal elevation with truncation. A higher grid reduces individual cell
size but does not make the threshold conservative.

The land model and atmosphere may limit what can be tiled immediately, so the
first deliverable is to measure lost/gained land and affected stores/fluxes on
every rung against the 10M support. The accepted architecture is either a
partial-cell/tile representation or an explicit model-form bracket that closes
all extensive terms. SPAT-5 owns that decision together with lake and
rootable/barren fractions.

### 4. Means are passed through nonlinear operators without the required moments

Rastetter et al.'s aggregation result applies directly. In general,
`process(mean(x))` is not `mean(process(x))`, and the error depends on curvature,
variance and covariance. The repository currently includes several such paths:

- `build_soil.py` mixes parent textures, then applies nonlinear weathering and
  pedogenesis at cell-mean climate;
- its `subgrid_slope` divides the within-cell elevation spread by native mesh
  spacing, so changing region count changes the inferred slope even when the
  generated physical relief has converged;
- `build_surface_roughness.py` averages a surface roughness and combines it with
  elevation spread before the model applies a logarithmic transfer coefficient;
  it also re-solves its orographic coefficient separately at each resolution to
  preserve a 2 m global mean, confounding a clean convergence comparison;
- the LPJ driver applies one forcing history and one soil column to all land in
  a climate cell; and
- feedback builders average cover or albedo before coupled snow, radiation,
  moisture and boundary-layer responses occur.

No generic sign follows for T42 versus T170. Productivity, evaporation or
temperature can move either way as variance, covariance, threshold crossings,
orographic precipitation and atmospheric feedback change. The old prediction
that finer resolution must lower total NPP from concavity is therefore
withdrawn. SPAT-4 defines semantic operators; SPAT-7 evaluates the material
nonlinear paths under both process-then-aggregate and aggregate-then-process.

### 5. One ecological column is not a spatial mosaic

Pedology retains continuous mixed texture in one soil map but the LPJ driver
chooses the dominant categorical soil code. Both collapse correlation among
lithology, elevation, soil depth, moisture access, lake/playa status and local
climate. A synthetic average soil can have hydraulic and nutrient behavior that
no part of the cell possesses, while a dominant class discards every minority
surface regardless of its ecological leverage.

LPJ's patches do not repair this. They sample stochastic demography under the
same cell environment; changing `npatch` or `patcharea` cannot create a wet
valley, exposed ridge or groundwater-fed fraction. SPAT-6 defines a limited set
of co-located ecological response units per climate cell with conserved area
weights, while DEMO-6 continues to own demographic sampling sensitivity.
LSHY-6, PLHY-5, BIO-11 and BIO-17 supply the hydrologic and surface fractions.

### 6. Finer climate support changes the whole coupled fixed point

The sequencing guide describes Loop D as a change after Loops A--C settle.
Although it then repeats some T85 work, the framing is wrong and the ceiling is
obsolete. Resolution changes terrain forcing and coastline, the atmospheric
solution and extremes, catchment precipitation/evaporation, recharge and water
table, pedogenesis and nutrients, vegetation, surface albedo/roughness and the
carve verdict. A T42 fixed point is not a fixed point of the T85--T170 system.

The detailed design in
`exoplasim/notes/restart-resolution-precision-converter.md` is the intended
bridge. It spectrally maps compatible atmospheric state, conservatively maps
gridpoint reservoirs and obtains target static fields from a target template.
Its output is explicitly a new initial condition, not bitwise continuation.
However, no open TASKS row owns the implementation despite the performance
roadmap naming it as a prerequisite. CLIM-49 closes that gap.

SPAT-8 then uses the converter for progressive advancement, but applies
convergence gates only after the candidate support has replayed every affected
coupled loop. If an intermediate rung is already converged for all declared
decisions, there is no obligation to run higher. If results have not converged
by T170, the result is a resolution bracket/structural uncertainty, not an
automatic endorsement of the highest answer.

### 7. Extensive biosphere products do not yet share a common area support

`score_prediction.py` reconstructs cell area and multiplies densities by the
full area of every binary land cell. It does not apply partial land or rootable
area. Similar ambiguity reaches soil-carbon and vegetation feedback products,
and a percentage of cells changes meaning when cell count and cell area change.

BIO-11 already owns the rootable fraction and BIO-12 the accepted equilibrium
window. SPAT-9 adds the spatial half: stock versus density semantics, native and
effective area, conservative remapping to one common comparison support, and
area-weighted pattern metrics. Resolution comparisons must not difference raw
arrays or totals defined over different coastlines.

### 8. No convergence protocol separates grid, aggregation and realization error

The existing Orogen resolution audit shows why a two-point comparison is not
enough: changing region count also perturbs the generator realization, and
several apparent terrain trends were smaller than same-resolution controls.
The climate ladder avoids that confound because all five exports share one
terrain hash, but it introduces others: changed land area, recalibrated
roughness, unequal equilibration and non-conservative outputs.

SPAT-8 pre-registers quantities and tolerances before execution. At minimum it
must assess water and energy closure; coastline and effective terrestrial area;
temperature/precipitation means, extremes and covariance; catchment water
balance and carve set; recharge/water-table classes; soil/hydraulic/nutrient
state; NPP, biomass, cover, biome area and carbon stocks; and surface feedback
fluxes. Every accepted candidate is compared on a common support after the same
equilibrium rule, with computational cost reported separately.

### 9. Support and operator provenance need fail-closed tests

Shape equality is not support equality. The repository has already experienced
half-world longitude shifts from matching coordinate labels across an
index-identical grid. A resolution ladder adds more opportunities to pair a
field, coupling matrix, restart, soil map or climatology with the wrong grid.

SPAT-1 makes support identity first-class; SPAT-10 enforces it. Tests should
cover constant preservation, area and mass closure, categorical fractions,
periodic longitude, Gaussian weights, target-mask changes, reduction identities,
operator order and round-trip/common-support comparison. Any nearest fallback,
discarded high mode, changed land/ocean ownership or unmapped extensive store
must appear in the report rather than be hidden by a finite output.

## Order of work

1. Land SPAT-1's support contract and SPAT-2's grid registry/path
   parameterization; add CLIM-49 as the owner of the already specified restart
   converter.
2. Extend GRID-2 through SPAT-3 to emit the 10M-under-grid support artifact for
   every requested rung, then implement SPAT-4's semantic/conservative
   operators and SPAT-10's no-simulation fixtures.
3. Close partial surface area under SPAT-5/BIO-11/BIO-17 and preserve the
   hydrologic/ecological mosaic under SPAT-6 with LSHY-6 and PLHY-5.
4. Audit and correct the material nonlinear operators under SPAT-7 rather than
   assuming a universal resolution-bias sign.
5. Make all feedback and scoring outputs support-aware under SPAT-9.
6. Only with explicit permission, use CLIM-49 to advance progressively through
   the T21/T42/T85/T127/T170 ladder and apply SPAT-8's pre-registered coupled
   convergence gates. Stop at the first sufficient support; carry a bracket if
   T170 is insufficient.
