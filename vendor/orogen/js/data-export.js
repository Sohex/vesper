// Raw data export — turns a generation result into serializable fields.
//
// DOM-free and Three.js-free so it runs identically in the browser and in Node
// (tools/export-planet.mjs). Two products:
//
//   raw   — one value per mesh region, in region order. Lossless: this is
//           exactly what the simulation computed, on the irregular Delaunay
//           sphere mesh it was computed on.
//   grid  — the same fields resampled onto a regular equirectangular lat/lon
//           grid, which is what ExoPlaSim / LPJ-GUESS and friends consume.
//
// COORDINATE CONVENTION. The whole app derives geography from the region
// position vector as lat = asin(y), lon = atan2(x, z) — note this is NOT the
// frame the Fibonacci point generator used internally (it built a z-up sphere).
// Everything user-visible (globe, map view, PNG exports, climate) uses the
// y-up convention, so the exporter does too. Longitude increases eastward from
// the +z axis; the grid runs -180..180 in lon and +90..-90 in lat (north-up,
// row 0 = north pole), matching the PNG exports.

import { KOPPEN_CLASSES as KOPPEN_TABLE } from './koppen.js';
import { BASIN_SURFACE_CLASSES, usesLandHeightBranch } from './basins.js';
import { ROCK_CLASSES, SHELF_SUBSTRATE_BOOTSTRAP } from './lithology.js';
import { elevToHeightKm, scaledHeightKm } from './color-map.js';
import { planetSummary } from './planet-params.js';
import { hashTypedArray, hashJson } from './sha256.js';
import {
    PLANET_RADIUS_KM, regionLatLon, regionCellArea,
    gaussianLatitudes, spectralGrid, latitudeEdges, rowForLatitude, gridCellArea,
} from './geometry.js';

export { PLANET_RADIUS_KM, regionLatLon, regionCellArea };

// ─────────────────────────────────────────────────────────────────────────
//  Geometry
// ─────────────────────────────────────────────────────────────────────────

/**
 * Per-region plate surface velocity, decomposed into east/north components.
 *
 * v = omega · (pole × p). omega is in the generator's arbitrary angular units,
 * not rad/Myr — the field is meaningful for direction and relative magnitude,
 * which is what the tectonic model itself uses it for.
 */
export function plateVelocityField(r_xyz, r_plate, plateVec) {
    const n = r_xyz.length / 3;
    const east = new Float32Array(n);
    const north = new Float32Array(n);
    const speed = new Float32Array(n);
    for (let r = 0; r < n; r++) {
        const pv = plateVec[r_plate[r]];
        if (!pv) continue;
        const x = r_xyz[3 * r], y = r_xyz[3 * r + 1], z = r_xyz[3 * r + 2];
        const [px, py, pz] = pv.pole;
        // v = omega · (pole × position)
        const vx = (py * z - pz * y) * pv.omega;
        const vy = (pz * x - px * z) * pv.omega;
        const vz = (px * y - py * x) * pv.omega;

        // Local east/north basis in the app's y-up frame
        let ex = z, ez = -x;
        const el = Math.sqrt(ex * ex + ez * ez);
        if (el > 1e-10) { ex /= el; ez /= el; } else { ex = 1; ez = 0; }
        let nx = y * ez, ny = z * ex - x * ez, nz = -y * ex;
        const nl = Math.sqrt(nx * nx + ny * ny + nz * nz) || 1;
        nx /= nl; ny /= nl; nz /= nl;

        east[r] = vx * ex + vz * ez;
        north[r] = vx * nx + vy * ny + vz * nz;
        speed[r] = Math.sqrt(vx * vx + vy * vy + vz * vz);
    }
    return { east, north, speed };
}

// ─────────────────────────────────────────────────────────────────────────
//  Field registry
//
//  Metadata for every exportable per-region field. Anything present in the
//  generation result but missing from this table is still exported, with a
//  generic descriptor — the point of this fork is that nothing stays internal.
// ─────────────────────────────────────────────────────────────────────────

const M = (units, description, extra = {}) => ({ units, description, ...extra });

export const FIELD_META = {
    // Geometry
    lat: M('degrees_north', 'Region centre latitude'),
    lon: M('degrees_east', 'Region centre longitude'),
    x: M('1', 'Region unit position vector, x component'),
    y: M('1', 'Region unit position vector, y component (sin latitude)'),
    z: M('1', 'Region unit position vector, z component'),
    cell_area: M('km2', 'Area of the region\'s dual mesh cell. A MESH property — raw output only, '
        + 'and the right weight for reducing MESH regions onto a grid cell. Resampled onto a grid it '
        + 'becomes the mean region area per cell, which carries no information about the cell, so it '
        + 'is not a grid weight; use grid_cell_area for that.'),
    grid_cell_area: M('km2', 'Area of this grid cell, exact, summing to the planet\'s surface area. '
        + 'The weight for area-averaging a gridded field. On a Gaussian grid the row boundaries are '
        + 'the Gauss-Legendre intervals, so this is grid/gauss_weights.bin times R² times the '
        + 'longitude step and a global mean taken with it is the one a spectral model takes; on an '
        + 'equally spaced grid the boundaries are midway in latitude.'),

    // Elevation
    // The model's internal elevation is a dimensionless shaping parameter, NOT
    // kilometres — it maps to height through a nonlinear hypsometric curve
    // (elevToHeightKm), where 0.5 is 1.1 km and 1.0 is 6 km. The parameter is
    // not bounded by 1 and the curve is linear at 6 km per unit above it, which
    // is where its shape function's domain ends. Anything consuming real
    // orography wants elevation_km.
    elevation: M('model units', 'Internal elevation parameter, NOT kilometres. Land maps to height '
        + 'through a nonlinear hypsometric curve; ocean is linear at 10 km per unit. Use '
        + 'elevation_km for physical height. Land is elevation > 0.'),
    elevation_km: M('km', 'Physical surface elevation relative to sea level, in kilometres. This is '
        + 'the orography a climate model wants. Land height includes the 1/g gravity scaling for '
        + 'this planet; ocean depth does not, because g cancels in the isostatic balance of the '
        + 'water column.'),
    elevation_pre_erosion: M('model units', 'Internal elevation parameter before terrain '
        + 'post-processing (erosion, detail noise, warp). Same nonlinear scale as elevation.'),
    land_mask: M('1', 'Elevation-sign test only: 1 where elevation > 0. This is NOT the land/sea '
        + 'mask — it floods every dry closed basin whose floor lies below sea level, which on this '
        + 'fork is the terrain the whole endorheic feature exists to preserve. Use surface_class '
        + 'for a land/sea mask; see manifest.landSeaMask for the size of the disagreement.',
        { categorical: true, dtype: 'uint8' }),

    // Plates
    plate_index: M('1', 'Dense plate index, 0..numPlates-1; see plates[] in the manifest', { categorical: true, dtype: 'int32' }),
    plate_seed_id: M('1', 'Raw plate ID as stored internally (a coarse-mesh seed region index)', { categorical: true, dtype: 'int32' }),
    plate_is_ocean: M('1', 'Whether the region\'s plate is oceanic', { categorical: true, dtype: 'uint8' }),
    super_plate_index: M('1', 'Super-plate (merged same-type plate group) index; -1 when super plates were not built', { categorical: true, dtype: 'int32' }),
    plate_velocity_east: M('model units', 'Eastward component of plate surface velocity from the Euler pole'),
    plate_velocity_north: M('model units', 'Northward component of plate surface velocity from the Euler pole'),
    plate_speed: M('model units', 'Magnitude of plate surface velocity'),

    // Tectonic state (stage 1)
    r_stress: M('1', 'Propagated tectonic stress magnitude'),
    r_stressDir_x: M('1', 'Stress direction unit vector, x'),
    r_stressDir_y: M('1', 'Stress direction unit vector, y'),
    r_stressDir_z: M('1', 'Stress direction unit vector, z'),
    r_subductFactor: M('1', 'Subduction propensity from plate density contrast; low values build mountains, high values build trenches'),
    r_boundaryType: M('1', 'Plate boundary class: 0 interior, 1 convergent, 2 divergent, 3 transform', { categorical: true, dtype: 'int8' }),
    r_bothOcean: M('1', 'Boundary cell where both plates are oceanic', { categorical: true, dtype: 'uint8' }),
    r_hasOcean: M('1', 'Boundary cell where at least one plate is oceanic', { categorical: true, dtype: 'uint8' }),
    r_mantleNorm: M('1', 'Mantle flow field, normalised to [-1, 1]'),
    mask_mountain: M('1', 'Seed cells for orogenic uplift', { categorical: true, dtype: 'uint8' }),
    mask_coastline: M('1', 'Seed cells for coastline elevation', { categorical: true, dtype: 'uint8' }),
    mask_oceanSeed: M('1', 'Seed cells for ocean floor elevation', { categorical: true, dtype: 'uint8' }),
    mask_stressMountain: M('1', 'Mountain seeds that collide rather than subduct', { categorical: true, dtype: 'uint8' }),

    // Spatial fields (stage 2) — hop counts, see bands{} for the resolution scaling
    r_isOcean: M('1', 'Ocean by plate type (distinct from land_mask, which is by elevation)', { categorical: true, dtype: 'uint8' }),
    dist_mountain: M('cells', 'BFS distance to nearest mountain seed'),
    dist_ocean: M('cells', 'BFS distance to nearest ocean seed'),
    dist_coastline: M('cells', 'BFS distance to nearest coastline seed'),
    dist_coast: M('cells', 'Ocean cells: BFS distance to the coast'),
    dist_coast_land: M('cells', 'Land cells: BFS distance to the coast through land'),
    dist_boundary: M('cells', 'BFS distance to the nearest plate boundary'),
    coastStressMax: M('1', 'Peak stress on the nearest stretch of coast'),
    coastSubductMax: M('1', 'Peak subduction factor on the nearest stretch of coast'),
    coastConvergent: M('1', 'Whether the nearest coast is on a convergent margin', { categorical: true, dtype: 'uint8' }),
    riftDist: M('cells', 'Distance into a continental rift zone; Infinity (raw) / NaN (grid) where no rift is in reach'),
    ridgeDist: M('cells', 'Distance from the mid-ocean ridge axis; Infinity (raw) / NaN (grid) where no ridge is in reach'),
    fractureDist: M('cells', 'Distance from the nearest oceanic fracture zone; Infinity (raw) / NaN (grid) where none is in reach'),
    backArcDist: M('cells', 'Distance into the back-arc basin behind a subduction margin'),
    backArcStress: M('1', 'Stress driving the back-arc extension'),

    // Terrain archetypes (stage 3)
    r_basinFactor: M('1', 'Sedimentary basin subsidence factor'),
    r_tectonicActivity: M('1', 'Composite tectonic activity used to modulate relief'),
    r_noiseAmp: M('1', 'Per-cell detail-noise amplitude multiplier'),
    r_t_foldBelt: M('1', 'Fold-belt archetype weight'),
    r_t_craton: M('1', 'Stable craton archetype weight'),
    r_t_basin: M('1', 'Basin archetype weight'),
    r_t_plateau: M('1', 'Plateau archetype weight'),

    // Elevation decomposition (debug layers)
    base: M('km', 'Base elevation contribution'),
    skeleton: M('km', 'Elevation after stage 4, before any noise or edifices — pure tectonic form'),
    tectonic: M('km', 'Tectonic uplift contribution'),
    noise: M('km', 'Fractal noise contribution'),
    noiseAmp: M('1', 'Noise amplitude field'),
    foldBeltWeight: M('1', 'Fold-belt weight (same as r_t_foldBelt)'),
    cratonWeight: M('1', 'Craton weight (same as r_t_craton)'),
    basinWeight: M('1', 'Basin weight (same as r_t_basin)'),
    basin: M('1', 'Basin subsidence factor'),
    interior: M('km', 'Plate-interior elevation contribution'),
    coastal: M('km', 'Coastal detail contribution'),
    ocean: M('km', 'Ocean floor contribution'),
    hotspot: M('km', 'Hotspot volcanism contribution'),
    lip: M('km', 'Large igneous province / flood basalt contribution'),
    tecActivity: M('1', 'Tectonic activity field'),
    margins: M('km', 'Continental margin contribution'),
    backArc: M('km', 'Back-arc basin contribution'),
    phasorRidge: M('km', 'Phasor-noise fold ridge contribution'),
    orogenicPower: M('1', 'Orogenic power, stored as (power - 0.5) in [-0.5, 0.5]'),
    uniformNoise: M('km', 'Uniform background land noise contribution'),
    dynamicTopo: M('km', 'Mantle-driven dynamic topography deflection'),
    erosionDelta: M('km', 'Net elevation change from the erosion / post-processing stage'),
    superPlates: M('1', 'Super-plate index as a float layer', { categorical: true }),

    // Plate physics diagnostics
    continentalDrag: M('1', 'Continental drag factor applied to the plate'),
    sizeVelocity: M('1', 'Size-velocity scaling factor applied to the plate'),
    plateSpeed: M('1', 'Plate angular velocity after all physics modifications'),
    velChange: M('1', 'Magnitude of velocity change caused by the physics pass'),
    mantleFlow: M('1', 'Mantle flow magnitude under the plate'),

    // Basins / hydrology
    basin_index: M('1', 'Index into basins.preserved[] of the preserved closed basin containing this region, or -1', { categorical: true, dtype: 'int32' }),
    is_endorheic: M('1', 'Region lies inside a preserved closed (endorheic) basin', { categorical: true, dtype: 'uint8' }),
    surface_class: M('1', 'Surface type: 0 ocean, 1 land, 2 inland water. AUTHORITATIVE for a '
        + 'land/sea mask — ocean means connected to the world ocean, so a dry closed basin below '
        + 'sea level is correctly land. Build an ExoPlaSim land/sea mask from this, not from '
        + 'land_mask.', { categorical: true, dtype: 'uint8' }),
    inland_water_basin: M('1', 'Index of the preserved basin whose lake covers this region, or -1', { categorical: true, dtype: 'int32' }),
    inland_water_level: M('km', 'Elevation of the lake surface covering this region; NaN where there is no inland water'),
    elevation_pre_conditioning: M('model units', 'Internal elevation parameter before any drainage '
        + 'conditioning — the surface the depression catalogue was measured on. Same nonlinear '
        + 'scale as `elevation`, NOT kilometres.'),
    drainage_terminal: M('1', 'Region index where this cell\'s flow path ends: an open-ocean cell, or '
        + 'a preserved basin\'s sink ON THE FINAL TERRAIN (basins[].finalCatchment.rootSink, which '
        + 'equals finalPreserved.sink — NOT the catalogue sink, which erosion may have left above '
        + 'the new low point). -2 if it resolves to neither. Grouping cell_area by this field '
        + 'reproduces finalCatchment.areaKm2 exactly; basins.drainageConsistency publishes the check.',
        { categorical: true, dtype: 'int32' }),

    // Lithology
    // Named `substrate_class`, not `surface_rock`, and the distinction is the
    // point: this is the top of the cover/basement stack, which is consolidated
    // lithology plus closed-basin fill. It is NOT what is at the surface.
    // Surficial cover -- soil, duricrust, loess, deflation armour, dune sand --
    // is a product of surface process history and is derived downstream. GLiM's
    // largest land class on Earth is unconsolidated sediment at 24.6%, and this
    // field deliberately has no counterpart for most of it.
    substrate_class: M('1', 'Substrate class: consolidated lithology plus closed-basin fill, at the top of the cover/basement stack. NOT the surface -- surficial cover is derived downstream. See rockClasses[] in the manifest', { categorical: true, dtype: 'uint8' }),
    basement_rock: M('1', 'Crystalline / oceanic basement rock class beneath the cover', { categorical: true, dtype: 'uint8' }),
    cover_rock: M('1', 'Sedimentary or volcanic cover rock class, where cover survives', { categorical: true, dtype: 'uint8' }),
    cover_thickness: M('km', 'Remaining thickness of the cover layer after erosion; 0 means basement is exposed'),
    substrate_class_pre_erosion: M('1', 'Substrate class before erosion, for comparison against substrate_class', { categorical: true, dtype: 'uint8' }),
    erodibility: M('1', 'Relative stream-power multiplier from the exposed rock, mean-normalised to 1 over land'),
    rock_albedo: M('1', 'Bare-rock shortwave albedo from the exposed rock class. A surface boundary '
        + 'condition for the climate stage; vegetation, snow AND INLAND WATER override it where '
        + 'they exist — this field is the rock, so a cell under a lake (surface_class == 2, see '
        + 'inland_water_level) still reports its dry-floor rock albedo and needs a water value '
        + 'substituted. Orogen does not apply it to anything.'),
    flow_accumulation: M('km2', 'Upstream contributing area draining through this region, in real area '
        + 'units so it compares across resolutions and planets. Routed on the depression-filled '
        + 'surface and terminated at basin SINKS, so a preserved basin\'s sink carries the whole '
        + 'inflow its lake receives.'),
    drain_to: M('1', 'Region index this cell drains to — steepest descent on the depression-filled '
        + 'surface, which crosses filled flats where the raw gradient is zero. -1 at an ocean cell '
        + 'or a preserved basin sink. Same tree as drainage_terminal and flow_accumulation.',
        { categorical: true, dtype: 'int32' }),
    ocean_basin_index: M('1', 'Marginal-sea index for submerged cells not connected to the world ocean, '
        + 'or -1; see hydrology.oceanBasins in the manifest for each basin\'s sill depth',
        { categorical: true, dtype: 'int32' }),
    scarp_potential: M('1', 'Relative likelihood that a sub-grid escarpment belongs here, 0..1. '
        + 'Orogen cannot resolve a cliff at this cell size; this marks where the ingredients are '
        + '(cover/basement erodibility contrast, proximity to where that interface daylights, and '
        + 'enough relief, meaning how far the ground stands above its own regional surroundings '
        + 'between two declared lengths) for a higher-resolution pass or a renderer to act on. '
        + 'Cover is always the weaker layer here, so these are stripped-edge plateau margins '
        + 'rather than resistant-caprock cuestas.'),

    // Climate (present only when climate was computed)
    tempSummer: M('degC', 'Mean summer surface temperature'),
    tempWinter: M('degC', 'Mean winter surface temperature'),
    tempContinentality: M('1', 'Continentality index used by the temperature model'),
    continentality: M('1', 'Continentality index from the wind model'),
    precipSummer: M('mm', 'Summer precipitation'),
    precipWinter: M('mm', 'Winter precipitation'),
    rainShadowSummer: M('1', 'Summer orographic rain shadow factor'),
    rainShadowWinter: M('1', 'Winter orographic rain shadow factor'),
    pressureSummer: M('hPa deviation', 'Summer sea-level pressure deviation'),
    pressureWinter: M('hPa deviation', 'Winter sea-level pressure deviation'),
    windSpeedSummer: M('model units', 'Summer wind speed'),
    windSpeedWinter: M('model units', 'Winter wind speed'),
    koppen: M('1', 'Köppen-Geiger class index; see koppenClasses in the manifest', { categorical: true, dtype: 'uint8' }),
    r_wind_east_summer: M('model units', 'Summer wind, eastward component'),
    r_wind_north_summer: M('model units', 'Summer wind, northward component'),
    r_wind_east_winter: M('model units', 'Winter wind, eastward component'),
    r_wind_north_winter: M('model units', 'Winter wind, northward component'),
    r_ocean_current_east_summer: M('model units', 'Summer ocean current, eastward component'),
    r_ocean_current_north_summer: M('model units', 'Summer ocean current, northward component'),
    r_ocean_current_east_winter: M('model units', 'Winter ocean current, eastward component'),
    r_ocean_current_north_winter: M('model units', 'Winter ocean current, northward component'),
    r_ocean_speed_summer: M('model units', 'Summer ocean current speed'),
    r_ocean_speed_winter: M('model units', 'Winter ocean current speed'),
    r_ocean_warmth_summer: M('1', 'Summer advected ocean warmth anomaly'),
    r_ocean_warmth_winter: M('1', 'Winter advected ocean warmth anomaly'),
    r_precip_summer: M('mm', 'Summer precipitation'),
    r_precip_winter: M('mm', 'Winter precipitation'),
    r_temperature_summer: M('degC', 'Summer temperature'),
    r_temperature_winter: M('degC', 'Winter temperature'),
};

/** Groups control which prefix a field lands under, and what --only can select. */
export const FIELD_GROUPS = {
    geometry: ['lat', 'lon', 'x', 'y', 'z', 'cell_area', 'grid_cell_area'],
    elevation: ['elevation', 'elevation_km', 'elevation_pre_erosion', 'land_mask'],
    plates: ['plate_index', 'plate_seed_id', 'plate_is_ocean', 'super_plate_index',
             'plate_velocity_east', 'plate_velocity_north', 'plate_speed',
             'continentalDrag', 'sizeVelocity', 'plateSpeed', 'velChange', 'mantleFlow'],
    tectonics: ['r_stress', 'r_stressDir_x', 'r_stressDir_y', 'r_stressDir_z', 'r_subductFactor',
                'r_boundaryType', 'r_bothOcean', 'r_hasOcean', 'r_mantleNorm',
                'mask_mountain', 'mask_coastline', 'mask_oceanSeed', 'mask_stressMountain',
                'r_isOcean', 'dist_mountain', 'dist_ocean', 'dist_coastline', 'dist_coast',
                'dist_coast_land', 'dist_boundary', 'coastStressMax', 'coastSubductMax',
                'coastConvergent', 'riftDist', 'ridgeDist', 'fractureDist', 'backArcDist', 'backArcStress',
                'r_basinFactor', 'r_tectonicActivity', 'r_noiseAmp',
                'r_t_foldBelt', 'r_t_craton', 'r_t_basin', 'r_t_plateau'],
    terrain: ['base', 'skeleton', 'tectonic', 'noise', 'noiseAmp', 'foldBeltWeight', 'cratonWeight',
              'basinWeight', 'basin', 'interior', 'coastal', 'ocean', 'hotspot', 'lip',
              'tecActivity', 'margins', 'backArc', 'phasorRidge', 'orogenicPower',
              'uniformNoise', 'dynamicTopo', 'erosionDelta', 'superPlates'],
    basins: ['basin_index', 'is_endorheic', 'surface_class', 'inland_water_basin',
             'inland_water_level', 'elevation_pre_conditioning', 'drainage_terminal'],
    orography: ['orog_mean', 'orog_std', 'orog_min', 'orog_max', 'orog_count',
                'orog_anisotropy', 'orog_angle'],
    lithology: ['substrate_class', 'basement_rock', 'cover_rock', 'cover_thickness',
                'substrate_class_pre_erosion', 'erodibility', 'scarp_potential', 'rock_albedo'],
    hydrology: ['flow_accumulation', 'drain_to', 'ocean_basin_index'],
    climate: ['tempSummer', 'tempWinter', 'tempContinentality', 'continentality',
              'precipSummer', 'precipWinter', 'rainShadowSummer', 'rainShadowWinter',
              'pressureSummer', 'pressureWinter', 'windSpeedSummer', 'windSpeedWinter', 'koppen',
              'r_wind_east_summer', 'r_wind_north_summer', 'r_wind_east_winter', 'r_wind_north_winter',
              'r_ocean_current_east_summer', 'r_ocean_current_north_summer',
              'r_ocean_current_east_winter', 'r_ocean_current_north_winter',
              'r_ocean_speed_summer', 'r_ocean_speed_winter',
              'r_ocean_warmth_summer', 'r_ocean_warmth_winter',
              'r_precip_summer', 'r_precip_winter', 'r_temperature_summer', 'r_temperature_winter'],
};

const GROUP_OF = (() => {
    const m = new Map();
    for (const [g, names] of Object.entries(FIELD_GROUPS)) for (const n of names) m.set(n, g);
    return m;
})();

function dtypeOf(arr) {
    if (arr instanceof Float32Array) return 'float32';
    if (arr instanceof Float64Array) return 'float64';
    if (arr instanceof Int32Array) return 'int32';
    if (arr instanceof Uint32Array) return 'uint32';
    if (arr instanceof Int16Array) return 'int16';
    if (arr instanceof Uint16Array) return 'uint16';
    if (arr instanceof Int8Array) return 'int8';
    if (arr instanceof Uint8Array) return 'uint8';
    return 'float32';
}

/**
 * Collect every per-region field from a generation result.
 *
 * @param {object} data  state.curData in the browser, or the pipeline ctx in Node
 * @returns {Array<{name, group, dtype, units, description, categorical, values}>}
 */
export function collectRegionFields(data) {
    const { mesh, r_xyz, t_xyz, r_elevation } = data;
    const n = mesh.numRegions;
    const fields = [];
    const seen = new Set();

    const add = (name, values, overrides = {}) => {
        if (!values || values.length !== n || seen.has(name)) return;
        seen.add(name);
        const meta = FIELD_META[name] || {};
        fields.push({
            name,
            group: overrides.group || GROUP_OF.get(name) || 'other',
            dtype: overrides.dtype || meta.dtype || dtypeOf(values),
            units: overrides.units ?? meta.units ?? 'unknown',
            description: overrides.description ?? meta.description ?? `Internal field "${name}" (undocumented)`,
            categorical: overrides.categorical ?? meta.categorical ?? false,
            values,
        });
    };

    // Geometry
    const { lat, lon } = regionLatLon(r_xyz);
    add('lat', lat);
    add('lon', lon);
    const xs = new Float32Array(n), ys = new Float32Array(n), zs = new Float32Array(n);
    for (let r = 0; r < n; r++) { xs[r] = r_xyz[3 * r]; ys[r] = r_xyz[3 * r + 1]; zs[r] = r_xyz[3 * r + 2]; }
    add('x', xs); add('y', ys); add('z', zs);
    // Prefer the pipeline's own cell areas: they were computed with this
    // planet's radius. Recomputing here without it silently returned Earth
    // areas, so a 1.2 R_earth planet reported Earth's surface area.
    if (data.cellArea) add('cell_area', data.cellArea);
    else if (t_xyz) add('cell_area', regionCellArea(mesh, r_xyz, t_xyz,
        data.planet ? data.planet.radiusKm : undefined));

    // Elevation. The model parameter and the physical height are different
    // things and both are exported; see FIELD_META for why.
    add('elevation', r_elevation);
    const reliefScale = data.planet ? data.planet.reliefScale : 1;
    // Land per surface_class, not per elevation sign. A dry closed-basin floor
    // below sea level is land, and running it through the bathymetric branch
    // read a -560 m salt pan as a -5.6 km trench.
    const surfaceCls = data.basins ? data.basins.r_surfaceClass : null;
    const elevKm = new Float32Array(n);
    for (let r = 0; r < n; r++) {
        const isLand = surfaceCls ? usesLandHeightBranch(surfaceCls[r]) : r_elevation[r] > 0;
        // Land scales by 1/g, ocean depth does not; scaledHeightKm carries
        // the reasoning. Do not "fix" the asymmetry.
        elevKm[r] = scaledHeightKm(r_elevation[r], reliefScale, isLand);
    }
    add('elevation_km', elevKm);
    add('elevation_pre_erosion', data.prePostElev);
    const landMask = new Uint8Array(n);
    for (let r = 0; r < n; r++) landMask[r] = r_elevation[r] > 0 ? 1 : 0;
    add('land_mask', landMask);

    // Plates
    if (data.r_plate) {
        add('plate_seed_id', data.r_plate instanceof Int32Array ? data.r_plate : Int32Array.from(data.r_plate));
        if (data.r_plateIndex) {
            add('plate_index', data.r_plateIndex);
        } else if (data.plateTable) {
            const idx = new Int32Array(n);
            for (let r = 0; r < n; r++) {
                const i = data.plateTable.seedToIndex.get(data.r_plate[r]);
                idx[r] = i === undefined ? -1 : i;
            }
            add('plate_index', idx);
        }
        if (data.plateIsOcean) {
            const isOcean = data.plateIsOcean instanceof Set ? data.plateIsOcean : new Set(data.plateIsOcean);
            const po = new Uint8Array(n);
            for (let r = 0; r < n; r++) po[r] = isOcean.has(data.r_plate[r]) ? 1 : 0;
            add('plate_is_ocean', po);
        }
        if (data.plateVec) {
            const v = plateVelocityField(r_xyz, data.r_plate, data.plateVec);
            add('plate_velocity_east', v.east);
            add('plate_velocity_north', v.north);
            add('plate_speed', v.speed);
        }
    }
    const sp = data.superPlates || data.superPlateData;
    if (sp && sp.r_superPlate) {
        add('super_plate_index', Int32Array.from(sp.r_superPlate));
    }

    // Tectonic internals
    const tect = data.tectonics;
    if (tect) {
        for (const [key, val] of Object.entries(tect)) {
            if (!val || typeof val !== 'object' || val.length === undefined) continue;
            if (key === 'r_stressDir') {
                const sx = new Float32Array(n), sy = new Float32Array(n), sz = new Float32Array(n);
                for (let r = 0; r < n; r++) { sx[r] = val[3 * r]; sy[r] = val[3 * r + 1]; sz[r] = val[3 * r + 2]; }
                add('r_stressDir_x', sx); add('r_stressDir_y', sy); add('r_stressDir_z', sz);
                continue;
            }
            add(key, val, { group: 'tectonics' });
        }
    }

    // Elevation decomposition + plate-physics diagnostics + climate rasters.
    // A few debug layers exist purely so the renderer can colour-map a field it
    // already has in a better-typed form; exporting those would ship the same
    // data twice, once as a lossy float copy.
    const RENDERER_DUPLICATES = new Set(['basinIndex', 'endorheic', 'surfaceClass',
        'surfaceRock', 'erodibility', 'coverThickness', 'scarpPotential',
        'flowAccum', 'oceanBasin']);
    if (data.debugLayers) {
        for (const [key, val] of Object.entries(data.debugLayers)) {
            if (!val || val.length !== n) continue;
            if (RENDERER_DUPLICATES.has(key)) continue;
            add(key, val);
        }
    }

    // Basins. These come from the pipeline's basin stage; absent when basin
    // preservation was disabled.
    if (data.basins) {
        const b = data.basins;
        add('basin_index', b.r_basinIndex);
        add('is_endorheic', b.r_isEndorheic);
        add('surface_class', b.r_surfaceClass);
        add('inland_water_basin', b.r_inlandWaterBasin);
        add('inland_water_level', b.r_inlandWaterLevelKm);
        add('elevation_pre_conditioning', b.preConditioningElev);
        if (b.r_finalTerminal) add('drainage_terminal', Int32Array.from(b.r_finalTerminal));
    }

    // Lithology. Absent when lithology was disabled.
    if (data.lithology) {
        const L = data.lithology;
        add('substrate_class', L.r_surfaceRock);
        add('basement_rock', L.r_basementRock);
        add('cover_rock', L.r_coverRock);
        add('cover_thickness', L.r_coverThicknessKm);
        add('substrate_class_pre_erosion', L.r_surfaceRockPreErosion);
        add('erodibility', L.r_erodibility);
        add('scarp_potential', L.r_scarpPotential);
        add('rock_albedo', L.r_albedo);
    }

    if (data.hydrology) {
        const H = data.hydrology;
        add('flow_accumulation', H.r_flowAccumKm2);
        add('drain_to', H.r_drainTo);
        add('ocean_basin_index', H.r_oceanBasinIndex);
    }

    // Climate vectors that live at the top level rather than in debugLayers
    for (const key of FIELD_GROUPS.climate) {
        if (key.startsWith('r_') && data[key]) add(key, data[key]);
    }

    return fields;
}

// ─────────────────────────────────────────────────────────────────────────
//  Equirectangular resampling
// ─────────────────────────────────────────────────────────────────────────

/**
 * Build a grid definition.
 *
 * type 'uniform' is equirectangular with equally spaced latitudes. type
 * 'gaussian' puts the rows at Gauss–Legendre latitudes, which is what spectral
 * models such as ExoPlaSim actually run on — emitting it directly from the mesh
 * removes a lossy equirect intermediate from the handoff. Longitudes are
 * equally spaced either way.
 */
export function makeGrid(width, height, type = 'uniform') {
    const lon = new Float64Array(width);
    for (let i = 0; i < width; i++) lon[i] = -180 + (i + 0.5) / width * 360;

    let lat, weights = null;
    if (type === 'gaussian') {
        const g = gaussianLatitudes(height);
        lat = g.lat;
        weights = g.weights;
    } else {
        lat = new Float64Array(height);
        for (let j = 0; j < height; j++) lat[j] = 90 - (j + 0.5) / height * 180;
    }
    return { width, height, type, lat, lon, weights, sinEdges: latitudeEdges(lat) };
}

/** Back-compat helper: cell centre coordinates for a uniform grid. */
export function gridCoords(width, height, type = 'uniform') {
    const g = makeGrid(width, height, type);
    return { lat: g.lat, lon: g.lon };
}

/** Row/column of the grid cell containing a given lat/lon. */
function cellOf(grid, latDeg, lonDeg) {
    const row = rowForLatitude(Math.sin(latDeg * Math.PI / 180), grid.sinEdges);
    let col = Math.floor((lonDeg + 180) / 360 * grid.width);
    if (col >= grid.width) col = grid.width - 1;
    if (col < 0) col = 0;
    return row * grid.width + col;
}

/**
 * Spatial index over region centres, bucketed by lat/lon, for nearest-region
 * queries. Bucket count targets ~2 regions per bucket so ring searches stay short.
 */
function buildRegionIndex(lat, lon) {
    const n = lat.length;
    const nLat = Math.max(2, Math.round(Math.sqrt(n / 4)));
    const nLon = nLat * 2;
    const counts = new Int32Array(nLat * nLon);
    const binOf = (r) => {
        let bi = Math.floor((90 - lat[r]) / 180 * nLat);
        if (bi >= nLat) bi = nLat - 1; if (bi < 0) bi = 0;
        let bj = Math.floor((lon[r] + 180) / 360 * nLon);
        if (bj >= nLon) bj = nLon - 1; if (bj < 0) bj = 0;
        return bi * nLon + bj;
    };
    for (let r = 0; r < n; r++) counts[binOf(r)]++;
    const offset = new Int32Array(nLat * nLon + 1);
    for (let b = 0; b < nLat * nLon; b++) offset[b + 1] = offset[b] + counts[b];
    const items = new Int32Array(n);
    const cursor = offset.slice(0, nLat * nLon);
    for (let r = 0; r < n; r++) items[cursor[binOf(r)]++] = r;
    return { nLat, nLon, offset, items };
}

/**
 * Nearest region for every cell of a width×height equirectangular grid.
 * Returns an Int32Array of region indices, row-major from the north pole.
 *
 * Nearest-centre lookup is the exact answer to "which cell does this point
 * fall in" for a Voronoi tessellation, so this is the correct sampler for
 * categorical fields, not just a convenient one.
 */
export function buildGridLookup(mesh, r_xyz, grid) {
    const { width, height } = grid;
    const { lat, lon } = regionLatLon(r_xyz);
    const idx = buildRegionIndex(lat, lon);
    const { nLat, nLon, offset, items } = idx;
    const out = new Int32Array(width * height);
    const n = lat.length;

    // Unit vectors for exact great-circle comparison (chord distance ordering
    // matches angular distance ordering, so no trig in the inner loop).
    const px = new Float32Array(n), py = new Float32Array(n), pz = new Float32Array(n);
    for (let r = 0; r < n; r++) { px[r] = r_xyz[3 * r]; py[r] = r_xyz[3 * r + 1]; pz[r] = r_xyz[3 * r + 2]; }

    const DEG = Math.PI / 180;
    for (let j = 0; j < height; j++) {
        const cellLat = grid.lat[j];
        const bi0 = Math.min(nLat - 1, Math.max(0, Math.floor((90 - cellLat) / 180 * nLat)));
        const clat = Math.cos(cellLat * DEG), slat = Math.sin(cellLat * DEG);
        for (let i = 0; i < width; i++) {
            const cellLon = grid.lon[i];
            const bj0 = Math.min(nLon - 1, Math.max(0, Math.floor((cellLon + 180) / 360 * nLon)));
            // Target point in the app's y-up frame: lat = asin(y), lon = atan2(x, z)
            const tx = clat * Math.sin(cellLon * DEG);
            const ty = slat;
            const tz = clat * Math.cos(cellLon * DEG);

            let best = -1, bestD = Infinity;
            // Expand the search ring until the ring's own minimum possible
            // distance exceeds the best hit found so far.
            for (let ring = 0; ring < nLat; ring++) {
                if (best >= 0) {
                    // A ring `ring` away is at least (ring-1) latitude bands off.
                    const minLatOff = (ring - 1) * (180 / nLat) * DEG;
                    if (minLatOff > 0 && 2 * Math.sin(minLatOff / 2) > Math.sqrt(bestD)) break;
                }
                let found = false;
                for (let di = -ring; di <= ring; di++) {
                    const bi = bi0 + di;
                    if (bi < 0 || bi >= nLat) continue;
                    const edge = (Math.abs(di) === ring);
                    for (let dj = -ring; dj <= ring; dj++) {
                        if (!edge && Math.abs(dj) !== ring) continue;
                        let bj = (bj0 + dj) % nLon; if (bj < 0) bj += nLon;
                        const b = bi * nLon + bj;
                        for (let k = offset[b], kEnd = offset[b + 1]; k < kEnd; k++) {
                            const r = items[k];
                            const dx = px[r] - tx, dy = py[r] - ty, dz = pz[r] - tz;
                            const d = dx * dx + dy * dy + dz * dz;
                            if (d < bestD) { bestD = d; best = r; }
                            found = true;
                        }
                    }
                }
                if (!found && ring > nLat) break;
            }
            out[j * width + i] = best;
        }
    }
    return out;
}

/**
 * Resample one region field onto an equirectangular grid.
 *
 * method 'nearest' takes the value of the region containing each grid-cell
 * centre. method 'mean' area-averages every region whose centre falls inside
 * the grid cell, falling back to nearest for cells no region centre lands in.
 * Downsampling (a coarse ExoPlaSim grid from a fine mesh) should use 'mean' for
 * continuous fields — nearest-sampling a 2.5M-cell mesh onto a 64×32 grid
 * would throw away 99.9% of the topography instead of averaging it.
 * Categorical fields always use nearest; averaging a class index is meaningless.
 */
export function resampleField(field, lookup, grid, lat, lon, cellArea, method) {
    const { width, height } = grid;
    const useMean = method === 'mean' && !field.categorical;
    const isFloat = field.dtype.startsWith('float') || useMean;
    const out = isFloat ? new Float32Array(width * height)
                        : makeTyped(field.dtype, width * height);

    if (!useMean) {
        for (let c = 0; c < out.length; c++) {
            const r = lookup[c];
            out[c] = r >= 0 ? field.values[r] : 0;
        }
        return out;
    }

    const sum = new Float64Array(width * height);
    const wsum = new Float64Array(width * height);
    const n = lat.length;
    let skippedNonFinite = 0;
    for (let r = 0; r < n; r++) {
        const v = field.values[r];
        // Several distance fields use Infinity for "no such feature within
        // reach". Averaging that would turn every grid cell containing one
        // region into Infinity, silently destroying the field — so non-finite
        // values are treated as missing and the cell averages what is left.
        if (!Number.isFinite(v)) { skippedNonFinite++; continue; }
        // Bin by the grid's own row edges, so this is correct for Gaussian
        // latitudes as well as equally spaced ones.
        const c = cellOf(grid, lat[r], lon[r]);
        const w = cellArea ? cellArea[r] : 1;
        sum[c] += v * w;
        wsum[c] += w;
    }
    for (let c = 0; c < out.length; c++) {
        if (wsum[c] > 0) { out[c] = sum[c] / wsum[c]; continue; }
        // Nothing finite fell in this cell: fall back to the containing region,
        // and emit NaN when that is non-finite too. NaN is the conventional
        // missing marker and survives numpy/xarray; Infinity does not.
        const r = lookup[c];
        const fallback = r >= 0 ? field.values[r] : 0;
        out[c] = Number.isFinite(fallback) ? fallback : NaN;
    }
    if (skippedNonFinite > 0) field.nonFiniteCount = skippedNonFinite;
    return out;
}

function makeTyped(dtype, len) {
    switch (dtype) {
        case 'float64': return new Float64Array(len);
        case 'int32': return new Int32Array(len);
        case 'uint32': return new Uint32Array(len);
        case 'int16': return new Int16Array(len);
        case 'uint16': return new Uint16Array(len);
        case 'int8': return new Int8Array(len);
        case 'uint8': return new Uint8Array(len);
        default: return new Float32Array(len);
    }
}

// ─────────────────────────────────────────────────────────────────────────
//  Sub-grid orography
// ─────────────────────────────────────────────────────────────────────────

export const SUBGRID_META = {
    orog_mean:   { units: 'km', description: 'Area-weighted mean elevation of the regions in this cell' },
    orog_std:    { units: 'km', description: 'Standard deviation of elevation within the cell — the standard sub-grid roughness input for surface drag and orographic gravity-wave schemes' },
    orog_min:    { units: 'km', description: 'Minimum region elevation within the cell' },
    orog_max:    { units: 'km', description: 'Maximum region elevation within the cell' },
    orog_count:  { units: '1',  description: 'Number of mesh regions whose centre falls in this cell; 0 means the cell was filled by nearest-region fallback and the statistics are undefined' },
    orog_anisotropy: { units: '1', description: 'Sub-grid slope anisotropy, 0 = perfectly ridged (all slopes share one axis) to 1 = isotropic. From the eigenvalues of the slope covariance tensor (Lott & Miller style)' },
    orog_angle:  { units: 'degrees', description: 'Orientation of the principal slope axis, measured anticlockwise from east; meaningful only where anisotropy is well below 1' },
};

/**
 * Statistics describing what the terrain does BELOW the grid.
 *
 * Resampling gives a cell its mean elevation and nothing else, but a climate
 * model wants to know how rough that cell is and whether the roughness is
 * organised into ridges. Both come out of the regions already assigned to each
 * cell, so this costs one extra pass.
 *
 * Slope is computed per region as a neighbour-averaged gradient in the local
 * east/north frame, then aggregated into a covariance tensor per cell. The
 * eigenvalues give anisotropy, the eigenvector gives the ridge orientation.
 *
 * `r_elevationKm` must be PHYSICAL height in km — the same `elevation_km` the
 * export publishes, hypsometric curve and 1/g relief scaling included. It used
 * to be handed the raw model parameter while SUBGRID_META declared km, which
 * made orog_mean bit-identical to the dimensionless `elevation` field and read
 * a 4.56 km peak as 1.14 "km". The error is not a constant factor, because the
 * model→km curve is nonlinear, so it could not be corrected downstream either.
 * orog_std is a drag and gravity-wave-scheme input, so this is a physical
 * quantity a climate model integrates, not a label.
 */
export function computeSubgridOrography(mesh, r_xyz, r_elevationKm, grid, lookup,
                                        lat, lon, cellArea, neighborDist, radiusKm) {
    const r_elevation = r_elevationKm;
    const { width, height } = grid;
    const nCells = width * height;
    const { adjOffset, adjList } = mesh;
    const n = mesh.numRegions;

    const sumW = new Float64Array(nCells);
    const sumZ = new Float64Array(nCells);
    const sumZ2 = new Float64Array(nCells);
    const minZ = new Float64Array(nCells).fill(Infinity);
    const maxZ = new Float64Array(nCells).fill(-Infinity);
    const count = new Float64Array(nCells);
    const sxx = new Float64Array(nCells);
    const syy = new Float64Array(nCells);
    const sxy = new Float64Array(nCells);

    const DEG = Math.PI / 180;
    for (let r = 0; r < n; r++) {
        const row = rowForLatitude(Math.sin(lat[r] * DEG), grid.sinEdges);
        let col = Math.floor((lon[r] + 180) / 360 * width);
        if (col >= width) col = width - 1; if (col < 0) col = 0;
        const c = row * width + col;

        const z = r_elevation[r];
        const w = cellArea ? cellArea[r] : 1;
        sumW[c] += w; sumZ[c] += z * w; sumZ2[c] += z * z * w;
        if (z < minZ[c]) minZ[c] = z;
        if (z > maxZ[c]) maxZ[c] = z;
        count[c]++;

        // Local east/north frame in the app's y-up convention.
        const x = r_xyz[3 * r], y = r_xyz[3 * r + 1], zz = r_xyz[3 * r + 2];
        let ex = zz, ez = -x;
        const el = Math.hypot(ex, ez);
        if (el > 1e-10) { ex /= el; ez /= el; } else { ex = 1; ez = 0; }
        let nx = y * ez, ny = zz * ex - x * ez, nz = -y * ex;
        const nl = Math.hypot(nx, ny, nz) || 1;
        nx /= nl; ny /= nl; nz /= nl;

        // Neighbour-averaged gradient, in rise over run (both in km).
        let gx = 0, gy = 0, deg = 0;
        for (let i = adjOffset[r], iEnd = adjOffset[r + 1]; i < iEnd; i++) {
            const nb = adjList[i];
            const dKm = (neighborDist ? neighborDist[i] : 0) * radiusKm;
            if (!(dKm > 0)) continue;
            const dz = (r_elevation[nb] - z) / dKm;
            let dx = r_xyz[3 * nb] - x, dy = r_xyz[3 * nb + 1] - y, dzv = r_xyz[3 * nb + 2] - zz;
            const dl = Math.hypot(dx, dy, dzv) || 1;
            dx /= dl; dy /= dl; dzv /= dl;
            gx += dz * (dx * ex + dzv * ez);
            gy += dz * (dx * nx + dy * ny + dzv * nz);
            deg++;
        }
        if (deg > 0) { gx /= deg; gy /= deg; }
        sxx[c] += gx * gx; syy[c] += gy * gy; sxy[c] += gx * gy;
    }

    const mean = new Float32Array(nCells);
    const std = new Float32Array(nCells);
    const mn = new Float32Array(nCells);
    const mx = new Float32Array(nCells);
    const cnt = new Float32Array(nCells);
    const aniso = new Float32Array(nCells);
    const angle = new Float32Array(nCells);

    let emptyCells = 0;
    for (let c = 0; c < nCells; c++) {
        cnt[c] = count[c];
        if (count[c] === 0) {
            // No region centre landed here. This is not missing data: the grid
            // cell is smaller than a mesh cell, so there is genuinely no
            // sub-grid variation to resolve — zero roughness is the correct
            // statement, and isotropic is the correct default orientation.
            // Emitting NaN instead would propagate through any drag or
            // gravity-wave scheme that consumes it.
            //
            // Unavoidable at high truncations near the poles, where Gaussian
            // rows are narrow while longitude spacing stays uniform: at T85 the
            // polar cells are 787 km2 against a 306 km2 mesh mean even at 2.5M
            // regions, so resolution does not remove them. orog_count stays 0 so
            // a consumer can still tell these cells apart.
            emptyCells++;
            const r = lookup[c];
            const h = r >= 0 ? r_elevation[r] : 0;
            mean[c] = h; mn[c] = h; mx[c] = h;
            std[c] = 0;
            aniso[c] = 1; angle[c] = 0;
            continue;
        }
        const m = sumZ[c] / sumW[c];
        mean[c] = m;
        const variance = Math.max(0, sumZ2[c] / sumW[c] - m * m);
        std[c] = Math.sqrt(variance);
        mn[c] = minZ[c];
        mx[c] = maxZ[c];

        // Eigen-decomposition of the 2x2 symmetric slope covariance.
        const a = sxx[c] / count[c], b = syy[c] / count[c], d = sxy[c] / count[c];
        const tr = a + b, det = a * b - d * d;
        const disc = Math.max(0, tr * tr / 4 - det);
        const l1 = tr / 2 + Math.sqrt(disc);          // major
        const l2 = tr / 2 - Math.sqrt(disc);          // minor
        aniso[c] = l1 > 1e-20 ? Math.sqrt(Math.max(0, l2) / l1) : 1;
        angle[c] = 0.5 * Math.atan2(2 * d, a - b) * 180 / Math.PI;
    }

    return {
        fields: {
            orog_mean: mean, orog_std: std, orog_min: mn, orog_max: mx,
            orog_count: cnt, orog_anisotropy: aniso, orog_angle: angle,
        },
        summary: {
            note: 'What the terrain does below the grid. orog_std is the usual roughness input '
                + 'for surface drag and orographic gravity-wave parameterisations; anisotropy and '
                + 'angle describe whether that roughness is ridged and along which axis.',
            emptyCells,
            emptyCellNote: 'Cells no region centre fell into, i.e. grid cells smaller than a mesh '
                + 'cell. Field values come from the containing region; orog_std is 0 and '
                + 'orog_anisotropy is 1, because a sub-mesh-scale cell resolves no sub-grid '
                + 'variation — that is a statement, not missing data. orog_count is 0 so these '
                + 'cells remain identifiable. Concentrated in the polar rows at high truncations, '
                + 'where Gaussian rows narrow but longitude spacing does not.',
        },
    };
}

// ─────────────────────────────────────────────────────────────────────────
//  Serialization
// ─────────────────────────────────────────────────────────────────────────

/** Bytes of a typed array, as a Uint8Array view (little-endian on all targets we run on). */
function bytesOf(arr) {
    return new Uint8Array(arr.buffer, arr.byteOffset, arr.byteLength);
}

/**
 * Build the complete export bundle.
 *
 * @param {object} data      generation result
 * @param {object} opts      { raw, grid, gridWidth, gridHeight, gridMethod, groups, onProgress }
 * @returns {{files: Array<{name: string, bytes: Uint8Array}>, manifest: object}}
 */
export function buildExportBundle(data, opts = {}) {
    const {
        raw = true,
        grid = true,
        gridWidth = 512,
        gridHeight = 256,
        gridMethod = 'mean',
        gridType = 'uniform',          // 'uniform' or 'gaussian'
        gridTruncation = null,         // e.g. 'T42', when the grid came from one
        subgridOrography = true,
        groups = null,
        onProgress = () => {},
    } = opts;

    onProgress(0.02, 'Collecting fields');
    const allFields = collectRegionFields(data);
    let fields = allFields;
    if (groups && groups.length) fields = fields.filter(f => groups.includes(f.group));
    // Area weighting needs cell areas whether or not the caller asked for them
    // in the output, so take them from the unfiltered set.
    const areaField = allFields.find(f => f.name === 'cell_area');

    const mesh = data.mesh;
    const files = [];
    const manifest = {
        format: 'world-orogen-data-export',
        version: 1,
        generatedBy: 'World Orogen (worldbuilding fork)',
        planetRadiusKm: PLANET_RADIUS_KM,
        seed: data.seed ?? null,
        params: data.params || data._params || null,
        numRegions: mesh.numRegions,
        numTriangles: mesh.numTriangles,
        coordinateConvention: {
            note: 'lat = asin(y), lon = atan2(x, z); longitude increases eastward from +z',
            latRange: [-90, 90],
            lonRange: [-180, 180],
        },
        elevation: {
            modelField: 'elevation',
            physicalField: 'elevation_km',
            units: 'km (elevation_km) / dimensionless (elevation)',
            seaLevel: 0,
            note: 'The `elevation` field is the model\'s internal shaping parameter, not '
                + 'kilometres: land maps to height through a nonlinear hypsometric curve where '
                + '0.5 is ~1.1 km and 1.0 is 6 km, and ocean is linear at 10 km per unit. The '
                + 'parameter is not bounded by 1; above it the curve is linear at 6 km per unit, '
                + 'which is where its shape function\'s domain ends. Use '
                + '`elevation_km` for physical orography. The land/ocean split is elevation > 0 '
                + 'in either field.',
            gravityScaling: 'elevation_km applies this planet\'s 1/g relief scaling to land '
                + 'heights; ocean depth is unscaled because g cancels in the isostatic balance.',
        },
        tectonicBands: data.tectonics ? data.tectonics.bands : null,
        boundaryTypes: { 0: 'interior', 1: 'convergent', 2: 'divergent', 3: 'transform' },
        koppenClasses: KOPPEN_CLASSES,
        plates: data.plates || (data.plateTable ? data.plateTable.plates : null),
        superPlates: summariseSuperPlates(data),
        basins: summariseBasins(data),
        surfaceClasses: BASIN_SURFACE_CLASSES,
        landSeaMask: summariseLandSeaMask(allFields),
        lithology: summariseLithology(data),
        planet: data.planetSummary || (data.planet ? planetSummary(data.planet) : null),
        hydrology: data.hydrology ? {
            note: 'Flow accumulation and ocean topology on the finished terrain. Discharge areas '
                + 'are the freshwater flux weighting a coastal ocean model wants.',
            // EVERY coastal land cell discharges to the sea, so the raw count is
            // a coastline length, not a river count — naming it riverMouthCount
            // asserted something untrue (median discharge is one cell's worth of
            // area). The distribution below is given instead of one arbitrary
            // threshold: pick the cut that suits the model being fed.
            coastalOutletCount: data.hydrology.mouthCount,
            outletNote: 'coastalOutletCount is every land cell that discharges to the open ocean '
                + '— effectively the coastline in cells. Use dischargeDistribution to pick a '
                + 'threshold for what counts as a river.',
            dischargeDistribution: data.hydrology.dischargeDistribution || null,
            largestRiverMouths: data.hydrology.mouths.slice(0, 50),
            oceanBasins: data.hydrology.oceanBasins,
            sillNote: 'A marginal sea\'s sill depth sets its character: deep sills exchange freely '
                + 'with the world ocean, shallow ones restrict it into Mediterranean-type '
                + '(evaporitic) or Black Sea-type (stratified, anoxic at depth) basins.',
        } : null,
        climateComputed: !!(data.debugLayers && data.debugLayers.koppen),
        raw: null,
        grid: null,
    };

    if (raw) {
        onProgress(0.1, 'Writing raw region fields');
        const entries = [];
        for (const f of fields) {
            const path = `raw/${f.name}.bin`;
            files.push({ name: path, bytes: bytesOf(f.values) });
            entries.push({
                name: f.name, group: f.group, path,
                dtype: f.dtype, shape: [mesh.numRegions],
                units: f.units, description: f.description, categorical: f.categorical,
            });
        }
        // Mesh topology, so a consumer can rebuild the dual mesh / walk neighbours.
        files.push({ name: 'raw/mesh_triangles.bin', bytes: bytesOf(toI32(mesh.triangles)) });
        files.push({ name: 'raw/mesh_halfedges.bin', bytes: bytesOf(toI32(mesh.halfedges)) });
        files.push({ name: 'raw/mesh_adj_offset.bin', bytes: bytesOf(toI32(mesh.adjOffset)) });
        files.push({ name: 'raw/mesh_adj_list.bin', bytes: bytesOf(toI32(mesh.adjList)) });
        manifest.raw = {
            layout: 'One value per mesh region, in region index order.',
            byteOrder: 'little-endian',
            fields: entries,
            mesh: {
                triangles: { path: 'raw/mesh_triangles.bin', dtype: 'int32', shape: [mesh.numSides],
                             description: 'Delaunay triangle corner region indices, 3 per triangle' },
                halfedges: { path: 'raw/mesh_halfedges.bin', dtype: 'int32', shape: [mesh.numSides],
                             description: 'Opposite half-edge index, -1 if none' },
                adjOffset: { path: 'raw/mesh_adj_offset.bin', dtype: 'int32', shape: [mesh.numRegions + 1],
                             description: 'CSR offsets into adjList' },
                adjList: { path: 'raw/mesh_adj_list.bin', dtype: 'int32', shape: [mesh.adjList.length],
                           description: 'Neighbour region indices, CSR-packed' },
            },
        };
    }

    if (grid) {
        onProgress(0.35, `Building ${gridWidth}×${gridHeight} ${gridType} grid lookup`);
        const grid = makeGrid(gridWidth, gridHeight, gridType);
        const lookup = buildGridLookup(mesh, data.r_xyz, grid);
        const { lat, lon } = regionLatLon(data.r_xyz);
        const cellArea = areaField ? areaField.values : null;
        const coords = { lat: grid.lat, lon: grid.lon };

        const entries = [];
        let manifestSubgrid = null;
        let done = 0;

        // The true area of each grid cell. Emitted before the resampled fields
        // so it is the obvious thing to reach for, and cell_area (a mesh
        // property) is excluded from the grid entirely — averaged onto a grid it
        // becomes the mean mesh-region area per cell, roughly constant, and
        // using it as a weight silently biases every integral toward the poles.
        {
            const gca = gridCellArea(grid, data.planet ? data.planet.radiusKm : PLANET_RADIUS_KM);
            files.push({ name: 'grid/grid_cell_area.bin', bytes: bytesOf(gca) });
            entries.push({
                name: 'grid_cell_area', group: 'geometry', path: 'grid/grid_cell_area.bin',
                dtype: 'float32', shape: [gridHeight, gridWidth],
                units: FIELD_META.grid_cell_area.units,
                description: FIELD_META.grid_cell_area.description,
                categorical: false, method: 'exact cell boundaries: Gauss-Legendre intervals on a '
                    + 'Gaussian grid, midway in latitude on an equally spaced one',
            });
        }

        for (const f of fields) {
            // Mesh-space areas are meaningless once resampled; see above.
            if (f.name === 'cell_area') continue;
            const values = resampleField(f, lookup, grid, lat, lon, cellArea, gridMethod);
            const path = `grid/${f.name}.bin`;
            files.push({ name: path, bytes: bytesOf(values) });
            entries.push({
                name: f.name, group: f.group, path,
                dtype: f.categorical ? f.dtype : (gridMethod === 'mean' ? 'float32' : f.dtype),
                shape: [gridHeight, gridWidth],
                units: f.units, description: f.description, categorical: f.categorical,
                method: f.categorical ? 'nearest' : gridMethod,
            });
            done++;
            onProgress(0.35 + 0.6 * (done / fields.length), `Resampling ${f.name}`);
        }
        files.push({ name: 'grid/lat.bin', bytes: bytesOf(new Float64Array(coords.lat)) });
        files.push({ name: 'grid/lon.bin', bytes: bytesOf(new Float64Array(coords.lon)) });
        if (grid.weights) {
            files.push({ name: 'grid/gauss_weights.bin', bytes: bytesOf(new Float64Array(grid.weights)) });
        }

        // Sub-grid orography. Only the cell MEAN of elevation survives normal
        // resampling, but a GCM's roughness and orographic gravity-wave drag
        // schemes want what the terrain does BELOW the grid: how rough it is,
        // and whether that roughness is ridged or isotropic.
        // Take elevation_km from the UNFILTERED set, the same way cell_area is:
        // these statistics are declared in km, so they must be computed from the
        // km field even when `--only orography` drops it from the output. Feeding
        // the raw model parameter here is what made orog_* dimensionless numbers
        // wearing a km label.
        const elevKmField = allFields.find(f => f.name === 'elevation_km');
        if (subgridOrography && elevKmField) {
            const stats = computeSubgridOrography(
                mesh, data.r_xyz, elevKmField.values, grid, lookup, lat, lon, cellArea,
                data.neighborDist, data.planet ? data.planet.radiusKm : PLANET_RADIUS_KM);
            for (const [name, arr] of Object.entries(stats.fields)) {
                const path = `grid/${name}.bin`;
                files.push({ name: path, bytes: bytesOf(arr) });
                entries.push({
                    name, group: 'orography', path,
                    dtype: 'float32', shape: [gridHeight, gridWidth],
                    units: SUBGRID_META[name].units,
                    description: SUBGRID_META[name].description,
                    categorical: false, method: 'subgrid statistic',
                });
            }
            manifestSubgrid = stats.summary;
        }
        manifest.grid = {
            projection: 'equirectangular (plate carrée), cell-centre registered',
            width: gridWidth, height: gridHeight,
            rowOrder: 'north to south (row 0 = +90)',
            areaWeighting: 'Use grid_cell_area for area-weighted statistics; on a Gaussian grid '
                + 'grid/gauss_weights.bin is the same partition expressed as weights summing to 2. '
                + 'The per-region cell_area field is raw-only and is NOT a valid grid weight.',
            colOrder: 'west to east (col 0 = -180)',
            defaultMethod: gridMethod,
            gridType,
            truncation: gridTruncation,
            latitudeNote: gridType === 'gaussian'
                ? 'Rows are at Gauss-Legendre latitudes, matching a spectral transform grid. They '
                + 'are quadrature abscissae and not cell centres, so a row owns its quadrature '
                + 'interval rather than the band around its midpoint; grid/gauss_weights.bin holds '
                + 'those intervals as weights summing to 2, and grid_cell_area is the same partition '
                + 'in km².'
                : 'Rows are equally spaced in latitude and are cell centres.',
            subgridOrography: manifestSubgrid,
            byteOrder: 'little-endian',
            coords: {
                lat: { path: 'grid/lat.bin', dtype: 'float64', shape: [gridHeight] },
                lon: { path: 'grid/lon.bin', dtype: 'float64', shape: [gridWidth] },
            },
            fields: entries,
        };
    }

    manifest.hashes = reproducibilityHashes(data, manifest);

    onProgress(0.98, 'Writing manifest');
    files.push({
        name: 'manifest.json',
        bytes: new TextEncoder().encode(JSON.stringify(manifest, null, 2)),
    });
    files.push({ name: 'README.txt', bytes: new TextEncoder().encode(readmeText(manifest)) });

    onProgress(1, 'Done');
    return { files, manifest, fields };
}

function toI32(arr) {
    return arr instanceof Int32Array ? arr : Int32Array.from(arr);
}

/**
 * Which field is the land/sea mask, and how much the alternatives disagree.
 *
 * Two land definitions coexist and they are not interchangeable. `land_mask` is
 * a bare elevation-sign test; `surface_class` asks whether a cell is connected
 * to the world ocean. They differ by exactly the terrain this fork exists to
 * preserve — dry closed-basin floors below sea level — so a consumer reaching
 * for the obvious-sounding name floods them. Quantified per export so the trap
 * is visible rather than latent.
 */
function summariseLandSeaMask(fields) {
    const get = (n) => { const f = fields.find(x => x.name === n); return f ? f.values : null; };
    const lm = get('land_mask'), sc = get('surface_class'), area = get('cell_area'),
          elev = get('elevation_km'), endo = get('is_endorheic');
    if (!lm || !sc) return null;

    let total = 0, byMask = 0, byClass = 0, disagree = 0, disagreeArea = 0;
    let minElev = Infinity, inBasin = 0;
    for (let r = 0; r < lm.length; r++) {
        const w = area ? area[r] : 1;
        total += w;
        if (lm[r] === 1) byMask += w;
        if (sc[r] === 1) byClass += w;
        if (sc[r] === 1 && lm[r] === 0) {
            disagree++; disagreeArea += w;
            if (elev && elev[r] < minElev) minElev = elev[r];
            if (endo && endo[r]) inBasin++;
        }
    }
    return {
        authoritative: 'surface_class',
        note: 'Build a land/sea mask from surface_class (0 ocean, 1 land, 2 inland water). '
            + 'land_mask is an elevation-sign test and floods dry closed basins that lie below '
            + 'sea level but are not connected to the ocean.',
        deepestDisagreeingElevationNote: 'In physical km, using the land branch of the elevation '
            + 'conversion — these are dry basin floors, not seabed.',
        landFractionBySurfaceClass: total > 0 ? byClass / total : 0,
        landFractionByLandMask: total > 0 ? byMask / total : 0,
        disagreementFraction: total > 0 ? disagreeArea / total : 0,
        disagreeingRegions: disagree,
        deepestDisagreeingElevationKm: minElev === Infinity ? null : minElev,
        disagreeingInsidePreservedBasin: inBasin,
        disagreementNote: 'Regions that surface_class calls land and land_mask calls ocean: dry '
            + 'closed-basin floor below sea level. Most sit inside a preserved basin; the rest are '
            + 'smaller enclosed depressions that did not clear the basin selection thresholds but '
            + 'are equally not sea.',
    };
}

/**
 * Basin section of the manifest.
 *
 * Deliberately reports three separable things, because conflating them is the
 * usual way endorheic modelling goes wrong:
 *   catalogue  — every depression the pre-conditioning terrain had (geometry)
 *   preserved  — the subset whose rims the drainage conditioning respected
 *   waterLevels — lake surfaces, only ever what the caller supplied
 *
 * Member region lists are omitted here (they can be millions of indices);
 * raw/basin_index.bin carries the per-region membership instead.
 */
function summariseBasins(data) {
    const b = data.basins;
    if (!b) return null;

    // Unsuffixed keys are the model's dimensionless elevation parameter; keys
    // ending Km / Km3 are physical. They used to be the same numbers under the
    // Km names, which read a -560 m salt pan as -5.6 km.
    const strip = (x) => ({
        id: x.id,
        index: x.index,
        parentIndex: x.parentIndex,
        nestDepth: x.nestDepth,
        spillsInto: x.spillsInto ?? null,
        sink: x.sink,
        sinkElevation: x.sinkElevation,
        spillElevation: x.spillElevation,
        depth: x.depth,
        sinkElevationKm: x.sinkElevationKm,
        spillElevationKm: x.spillElevationKm,
        depthKm: x.depthKm,
        areaKm2: x.areaKm2,
        volumeKm3: x.volumeKm3,
        cellCount: x.cellCount,
        spillFrom: x.spillFrom,
        spillTo: x.spillTo,
    });

    const preserved = b.selected.map((x, i) => {
        const fin = b.finalMeasurements[i] || {};
        const cat = b.finalCatchments[i] || {};
        const drift = b.catchmentDrift[i] || {};
        const water = b.waterLevels[i] || {};
        // The two records must name the same cell. `finalPreserved.sink` comes
        // from remeasureBasins and `finalCatchment.rootSink` from the routing;
        // when they drifted apart, drainage_terminal and finalCatchment.areaKm2
        // shipped as two contradictory answers to the same question and nothing
        // on the artifact said so. Cheap to check, so it is checked here rather
        // than trusted.
        if (fin.sink !== undefined && cat.sink !== undefined && fin.sink !== cat.sink) {
            throw new Error(`Basin ${x.id}: finalPreserved.sink ${fin.sink} != catchment root `
                + `${cat.sink} — the drainage field and the per-basin area would disagree`);
        }
        return {
            ...strip(x),
            selectedBy: x.selectedBy,
            hypsometry: x.hypsometry,
            natural: {
                spillDepth: x.depth,
                spillDepthKm: x.depthKm,
                volumeKm3: x.volumeKm3,
                areaKm2: x.areaKm2,
            },
            finalPreserved: {
                sink: fin.sink,
                sinkElevation: fin.sinkElevation,
                spillElevation: fin.spillElevation,
                spillDepth: fin.spillDepth,
                sinkElevationKm: fin.sinkElevationKm,
                spillElevationKm: fin.spillElevationKm,
                spillDepthKm: fin.spillDepthKm,
                retainedFraction: fin.retainedFraction,
                volumeKm3: fin.volumeKm3,
                floodedAreaKm2: fin.floodedAreaKm2,
            },
            finalCatchment: {
                areaKm2: cat.catchmentAreaKm2,
                basinFloorAreaKm2: cat.basinFloorAreaKm2,
                catchmentToFloorRatio: cat.catchmentToFloorRatio,
                // The cell drainage_terminal points at for every cell in this
                // catchment. Equals finalPreserved.sink, and differs from the
                // catalogue sink whenever erosion moved the basin's low point.
                rootSink: cat.sink,
                catalogueSink: cat.catalogueSink,
                measuredOn: 'depression-filled routing — every land cell resolves to a terminal, '
                          + 'so this is the full contributing area to integrate precipitation over. '
                          + 'Group drainage_terminal by rootSink (= finalPreserved.sink) to '
                          + 'reproduce areaKm2 exactly; grouping by the catalogue sink will not, '
                          + 'because that cell is not a root on the final terrain.',
            },
            retain: x.retain ?? 1,
            catchmentDriftFromPreConditioning: {
                jaccard: drift.jaccard,
                addedCells: drift.addedCount,
                removedCells: drift.removedCount,
            },
            waterLevelKm: water.levelKm ?? null,
        };
    });

    return {
        units: 'Unsuffixed elevation/depth keys are the model\'s dimensionless elevation '
             + 'parameter. Keys ending Km or Km3 are physical, converted through the LAND branch '
             + 'of the height curve — a closed-basin floor below sea level is dry ground, not '
             + 'seabed. retainedFraction is a ratio of model-unit depths, so the nonlinear curve '
             + 'does not distort it.',
        note: 'Depression geometry is measured here; whether a basin is endorheic is a water '
            + 'balance and is not decided by Orogen. hypsometry gives level/area/volume for that '
            + 'calculation, finalCatchment gives the area to integrate precipitation over.',
        resolution: b.resolution,
        // Published proof that the gridded/raw field and the per-basin table
        // were derived from the same routing and the same sink definition.
        drainageConsistency: b.drainageConsistency ?? null,
        sharedCatchmentRoots: b.sharedCatchmentRoots ?? [],
        // Which drainage hypothesis produced this export. Two runs of the same
        // seed now differ by the carve decision alone, and the elevation hash
        // does not say why — this does.
        selectionSource: b.selectionSource ?? 'threshold',
        drainageHypothesis: b.drainageHypothesis ?? null,
        incision: b.incision ? { ...b.incision,
            note: 'Basins with retain < 1 had a notch cut at their saddle, lowering the rim by '
                + '(1 - retain) of the basin relief. Merely permitting the carve to go lower does '
                + 'nothing — with a terminal root the flood has no deficit there, so the cut is '
                + 'made deliberately.' } : null,
        catalogueStability: 'Basin ids are computed from the PRE-CONDITIONING surface, which is '
            + 'produced before any preserve/carve decision. They are therefore identical across '
            + 'carve iterations of the same seed, parameters and region count — which is what lets '
            + 'a water-balance loop send verdicts back by id. Verify with hashes.basinCatalogue: '
            + 'if it matches, a verdict list computed against an earlier export is still valid.',
        selectionCriteria: {
            minDepthKm: b.resolution.minDepthKm,
            minAreaKm2: b.resolution.minAreaKm2,
            minCells: b.resolution.minCells,
            rationale: 'Signal-vs-noise floors: is the depression a landform the mesh resolves? '
                     + 'Not a judgement about endorheism.',
        },
        counts: {
            depressionsDetected: b.catalogue.length,
            preserved: b.selected.length,
            belowSeaLevelComponentsKept: b.topologyFixup ? b.topologyFixup.keptComponents : null,
            belowSeaLevelComponentsFilled: b.topologyFixup ? b.topologyFixup.filledComponents : null,
        },
        warnings: b.warnings,
        preserved,
        catalogue: b.catalogue.map(strip),
    };
}

/**
 * Lithology section of the manifest.
 *
 * Parent material only — this is not soil. Soil needs climate, vegetation and
 * time, so it belongs to a stage downstream of ExoPlaSim that combines this
 * rock map with a real climate.
 */
function summariseLithology(data) {
    const L = data.lithology;
    if (!L) return null;
    return {
        note: 'Parent material, not soil. Rock class is derived from tectonic history alone; '
            + 'soil texture requires climate and belongs to a later stage. Erosion responded to '
            + 'these rock types via the erodibility field.',
        rockClasses: ROCK_CLASSES,
        erodibility: {
            strength: L.strength ?? null,
            meanBeforeNormalisation: L.erodibilityMeanBeforeNormalisation,
            note: 'Mean-normalised to 1 over land before erosion, so lithology redistributes '
                + 'erosion without changing its total. Values drift slightly from 1 afterwards '
                + 'as cells exhume to basement.',
        },
        exhumedCells: L.exhumedCells,
        scarp: {
            cellsAbove025: L.scarpCells,
            note: 'scarp_potential marks where an escarpment belongs, not geometry. A cliff is a '
                + 'sub-kilometre feature and cells here are tens of km across, so the landform '
                + 'lives below the grid entirely. Intended for a higher-resolution downstream '
                + 'pass or a renderer.',
            kind: 'Cover is always weaker than basement in this model, so these are stripped-edge '
                + 'scarps (plateau margins, shield/cover boundaries) rather than resistant-caprock '
                + 'cuestas, which would need a third layer.',
        },
        exhumationNote: 'Cells whose exposed rock changed during erosion because the cover layer '
                      + 'was stripped away, revealing basement.',
        // A LABELLED BOOTSTRAP, not a result. This is the one climate call in a
        // module whose own planetSummary says Orogen does not do climate, and it
        // is made with an Earth constant. Emitted as structure rather than as a
        // sentence so a downstream reclassifier can key on it, and declared
        // beside the rule in lithology.js so the label cannot drift from it.
        shelfSubstrateBootstrap: SHELF_SUBSTRATE_BOOTSTRAP,
        compositionLand: L.composition,
        compositionGlobal: L.compositionSeafloor,
        // State the denominator. compositionLand used to be measured against
        // elevation > 0, which excludes the dry sub-sea-level basin floors and
        // under-reported evaporite — the class that accumulates on exactly those
        // floors — by 2.2 points. Both fractions here are area-weighted over the
        // named set.
        compositionDenominator: {
            land: L.landDenominator ?? null,
            landAreaKm2: L.landAreaKm2 ?? null,
            global: 'every region, land and seafloor',
            note: 'This export has two land definitions — see manifest.landSeaMask. Any published '
                + 'fraction must name which one it used; these do.',
        },
    };
}

/**
 * Canonical content hashes, so a downstream run can assert it is looking at the
 * planet it thinks it is. Same seed and parameters ⇒ same hashes, in the
 * browser and headless alike.
 */
function reproducibilityHashes(data, manifest) {
    const h = {};
    if (data.r_elevation) h.finalElevation = hashTypedArray(data.r_elevation);
    if (data.basins && data.basins.preConditioningElev) {
        h.preConditioningElevation = hashTypedArray(data.basins.preConditioningElev);
    }
    if (data.prePostElev) h.preErosionElevation = hashTypedArray(data.prePostElev);
    if (manifest.basins) {
        h.basinCatalogue = hashJson(manifest.basins.catalogue);
        h.basinsPreserved = hashJson(manifest.basins.preserved.map(b => b.id));
    }
    if (data.basins && data.basins.r_basinIndex) {
        h.basinMembership = hashTypedArray(data.basins.r_basinIndex);
    }
    if (data.lithology) {
        h.surfaceRock = hashTypedArray(data.lithology.r_surfaceRock);
        h.erodibility = hashTypedArray(data.lithology.r_erodibility);
    }
    if (manifest.params) h.params = hashJson(manifest.params);
    return h;
}

function summariseSuperPlates(data) {
    const sp = data.superPlates || data.superPlateData;
    if (!sp) return null;
    const isOcean = sp.superPlateIsOcean instanceof Set
        ? sp.superPlateIsOcean : new Set(sp.superPlateIsOcean || []);
    const out = [];
    for (let i = 0; i < sp.numSuperPlates; i++) {
        const pv = sp.superPlateVec ? sp.superPlateVec[i] : null;
        out.push({
            index: i,
            isOcean: isOcean.has(i),
            density: sp.superPlateDensity ? sp.superPlateDensity[i] ?? null : null,
            eulerPole: pv ? [pv.pole[0], pv.pole[1], pv.pole[2]] : null,
            omega: pv ? pv.omega : null,
        });
    }
    return out;
}

/** Köppen class index → code, in the order classifyKoppen emits. */
export const KOPPEN_CLASSES = KOPPEN_TABLE.map((c, i) => ({ id: i, code: c.code, name: c.name }));

function readmeText(manifest) {
    return `World Orogen data export
========================

Seed: ${manifest.seed}
Regions: ${manifest.numRegions}
${manifest.grid ? `Grid: ${manifest.grid.width} x ${manifest.grid.height} equirectangular` : 'Grid: not exported'}

Everything is described in manifest.json: for each field you get its path,
dtype, shape, units and a description. All binaries are flat little-endian
arrays with no header.

Reading a field in Python:

    import json, numpy as np
    m = json.load(open('manifest.json'))
    f = next(x for x in m['grid']['fields'] if x['name'] == 'elevation')
    elev = np.fromfile(f['path'], dtype=f['dtype']).reshape(f['shape'])
    # elev[0] is the northernmost row, elev[:, 0] is longitude -180.

Coordinates: ${manifest.coordinateConvention.note}
Elevation: ${manifest.elevation.units}, sea level at ${manifest.elevation.seaLevel}.
${manifest.elevation.note}

raw/  — one value per mesh region (irregular Delaunay sphere mesh). Region
        positions are in raw/lat.bin, raw/lon.bin and raw/x|y|z.bin; cell areas
        in raw/cell_area.bin. Mesh topology is included so the dual mesh can be
        rebuilt and neighbours walked.
grid/ — the same fields resampled to a regular lat/lon grid. Continuous fields
        are area-weighted means of the regions falling in each cell; categorical
        fields (plate index, boundary type, Köppen class, masks) take the value
        of the region containing the cell centre.
`;
}

// ─────────────────────────────────────────────────────────────────────────
//  ZIP (stored, no compression) — used by the browser download path
// ─────────────────────────────────────────────────────────────────────────

const CRC_TABLE = (() => {
    const t = new Uint32Array(256);
    for (let i = 0; i < 256; i++) {
        let c = i;
        for (let k = 0; k < 8; k++) c = (c & 1) ? (0xEDB88320 ^ (c >>> 1)) : (c >>> 1);
        t[i] = c >>> 0;
    }
    return t;
})();

function crc32(bytes) {
    let c = 0xFFFFFFFF;
    for (let i = 0; i < bytes.length; i++) c = CRC_TABLE[(c ^ bytes[i]) & 0xFF] ^ (c >>> 8);
    return (c ^ 0xFFFFFFFF) >>> 0;
}

/**
 * Pack files into a ZIP archive with no compression.
 *
 * Zip64 is used unconditionally for the per-entry records so exports larger
 * than 4 GB (a 2.5M-region mesh with every field is easily that) don't silently
 * truncate. Stored entries keep the payload byte-identical to the .bin files.
 */
export function zipStore(files) {
    const enc = new TextEncoder();
    const locals = [];
    const central = [];
    let offset = 0;

    for (const f of files) {
        const nameBytes = enc.encode(f.name);
        const data = f.bytes;
        const crc = crc32(data);

        // Local file header + Zip64 extra field
        const lh = new Uint8Array(30 + nameBytes.length + 20);
        const lv = new DataView(lh.buffer);
        lv.setUint32(0, 0x04034b50, true);
        lv.setUint16(4, 45, true);           // version needed (4.5 = zip64)
        lv.setUint16(6, 0, true);            // flags
        lv.setUint16(8, 0, true);            // method = stored
        lv.setUint16(10, 0, true);           // time
        lv.setUint16(12, 0x21, true);        // date (1996-01-01, fixed for reproducibility)
        lv.setUint32(14, crc, true);
        lv.setUint32(18, 0xFFFFFFFF, true);  // compressed size in zip64 extra
        lv.setUint32(22, 0xFFFFFFFF, true);  // uncompressed size in zip64 extra
        lv.setUint16(26, nameBytes.length, true);
        lv.setUint16(28, 20, true);          // extra field length
        lh.set(nameBytes, 30);
        const ev = new DataView(lh.buffer, 30 + nameBytes.length);
        ev.setUint16(0, 0x0001, true);
        ev.setUint16(2, 16, true);
        ev.setBigUint64(4, BigInt(data.length), true);
        ev.setBigUint64(12, BigInt(data.length), true);

        locals.push(lh, data);
        central.push({ nameBytes, crc, size: data.length, offset });
        offset += lh.length + data.length;
    }

    const centralStart = offset;
    const cdParts = [];
    for (const c of central) {
        const ch = new Uint8Array(46 + c.nameBytes.length + 28);
        const cv = new DataView(ch.buffer);
        cv.setUint32(0, 0x02014b50, true);
        cv.setUint16(4, 45, true);
        cv.setUint16(6, 45, true);
        cv.setUint16(8, 0, true);
        cv.setUint16(10, 0, true);
        cv.setUint16(12, 0, true);
        cv.setUint16(14, 0x21, true);
        cv.setUint32(16, c.crc, true);
        cv.setUint32(20, 0xFFFFFFFF, true);
        cv.setUint32(24, 0xFFFFFFFF, true);
        cv.setUint16(28, c.nameBytes.length, true);
        cv.setUint16(30, 28, true);
        cv.setUint16(32, 0, true);
        cv.setUint16(34, 0, true);
        cv.setUint16(36, 0, true);
        cv.setUint32(38, 0, true);
        cv.setUint32(42, 0xFFFFFFFF, true);
        ch.set(c.nameBytes, 46);
        const ev = new DataView(ch.buffer, 46 + c.nameBytes.length);
        ev.setUint16(0, 0x0001, true);
        ev.setUint16(2, 24, true);
        ev.setBigUint64(4, BigInt(c.size), true);
        ev.setBigUint64(12, BigInt(c.size), true);
        ev.setBigUint64(20, BigInt(c.offset), true);
        cdParts.push(ch);
        offset += ch.length;
    }
    const centralSize = offset - centralStart;

    // Zip64 end of central directory + locator + classic EOCD
    const end = new Uint8Array(56 + 20 + 22);
    const dv = new DataView(end.buffer);
    dv.setUint32(0, 0x06064b50, true);
    dv.setBigUint64(4, BigInt(44), true);
    dv.setUint16(12, 45, true);
    dv.setUint16(14, 45, true);
    dv.setUint32(16, 0, true);
    dv.setUint32(20, 0, true);
    dv.setBigUint64(24, BigInt(central.length), true);
    dv.setBigUint64(32, BigInt(central.length), true);
    dv.setBigUint64(40, BigInt(centralSize), true);
    dv.setBigUint64(48, BigInt(centralStart), true);
    dv.setUint32(56, 0x07064b50, true);
    dv.setUint32(60, 0, true);
    dv.setBigUint64(64, BigInt(offset), true);
    dv.setUint32(72, 1, true);
    dv.setUint32(76, 0x06054b50, true);
    dv.setUint16(80, 0xFFFF, true);
    dv.setUint16(82, 0xFFFF, true);
    dv.setUint16(84, 0xFFFF, true);
    dv.setUint16(86, 0xFFFF, true);
    dv.setUint32(88, 0xFFFFFFFF, true);
    dv.setUint32(92, 0xFFFFFFFF, true);
    dv.setUint16(96, 0, true);

    const total = offset + end.length;
    const out = new Uint8Array(total);
    let p = 0;
    for (const part of locals) { out.set(part, p); p += part.length; }
    for (const part of cdParts) { out.set(part, p); p += part.length; }
    out.set(end, p);
    return out;
}
