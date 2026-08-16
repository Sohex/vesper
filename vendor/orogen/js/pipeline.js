// Headless generation pipeline.
//
// This is the single implementation of "build a planet from parameters". It
// has no DOM, no Worker globals and no CDN imports, so it runs unchanged in
// the browser worker (js/planet-worker.js) and in Node (tools/export-planet.mjs).
//
// The caller is responsible for calling setDelaunator() before invoking
// runGeneratePipeline — the browser worker loads Delaunator from a CDN, Node
// loads it from node_modules.

import { makeRng } from './rng.js';
import { SimplexNoise } from './simplex-noise.js';
import { buildSphere, generateTriangleCenters, computeNeighborDist } from './sphere-mesh.js';
import { generateCoarsePlates, projectCoarsePlates } from './coarse-plates.js';
import { smoothAndReconnectPlates } from './plates.js';
import { assignElevation } from './elevation.js';
import { buildSuperPlates } from './super-plates.js';
import { warpTerrain, smoothElevation, erodeComposite, sharpenRidges, applySoilCreep, applyDetailNoise } from './terrain-post.js';
import { computeWind } from './wind.js';
import { computeOceanCurrents } from './ocean.js';
import { computePrecipitation } from './precipitation.js';
import { computeTemperature } from './temperature.js';
import { classifyKoppen } from './koppen.js';
import { applyPlatePhysics, expandPlatePhysicsDebug } from './plate-physics.js';
import {
    detectBasins, selectBasins, buildBasinProtection, remeasureBasins,
    computeFinalPreservedCatchments, compareCatchments, applyInlandWaterLevels,
    surfaceClasses, basinResolutionContext, findOpenOcean, attachHypsometry, parseBasinList,
    inciseOutlets, auditCatchmentConsistency,
    computeRivers, computeOceanBasins,
} from './basins.js';
import { regionCellArea } from './geometry.js';
import { makePlanet, planetSummary } from './planet-params.js';
import {
    classifyLithology, buildLithoState, finalSurfaceRock, rockComposition, ROCK_CLASSES,
    computeScarpPotential, buildAlbedo, saltCrustMask,
} from './lithology.js';
import { SUPER_PLATE_PHYSICS_MULT, DETAIL_NOISE_DAMPEN_STRENGTH } from './terrain-config.js';

// Compute triangle elevations from region elevations
export function computeTriangleElevations(mesh, r_elevation) {
    const t_elevation = new Float32Array(mesh.numTriangles);
    for (let t = 0; t < mesh.numTriangles; t++) {
        const s0 = 3 * t;
        const a = mesh.s_begin_r(s0), b = mesh.s_begin_r(s0 + 1), c = mesh.s_begin_r(s0 + 2);
        t_elevation[t] = (r_elevation[a] + r_elevation[b] + r_elevation[c]) / 3;
    }
    return t_elevation;
}

// Combined craton/basin dampen field for detail noise (1 = max dampen).
// Returns null when geological annotations aren't available (e.g. heightmap imports).
export function computeDetailDampenField(debugLayers) {
    const cw = debugLayers && debugLayers.cratonWeight;
    const bw = debugLayers && debugLayers.basinWeight;
    if (!cw || !bw) return null;
    const N = cw.length;
    const r_dampen = new Float32Array(N);
    for (let r = 0; r < N; r++) {
        const a = cw[r], b = bw[r];
        r_dampen[r] = a > b ? a : b;
    }
    return r_dampen;
}

// Orogenic power as a [0, 1] amplitude multiplier for detail noise.
// debugLayers.orogenicPower is stored as raw_oroPower − 0.5 (for diverging
// colormap), so we add 0.5 and clamp to recover the [0, 1] factor.
export function computeOrogenicField(debugLayers) {
    const op = debugLayers && debugLayers.orogenicPower;
    if (!op) return null;
    const N = op.length;
    const r_oro = new Float32Array(N);
    for (let r = 0; r < N; r++) {
        const v = op[r] + 0.5;
        r_oro[r] = v < 0 ? 0 : (v > 1 ? 1 : v);
    }
    return r_oro;
}

// Run terrain post-processing with per-step timing
export function runPostProcessing(mesh, r_xyz, r_elevation, params, neighborDist, seed, r_hotspot, r_dampen, r_orogenic, basinOpts = null, lithoOpts = null, hydro = null) {
    const { smoothing, glacialErosion, hydraulicErosion, thermalErosion, ridgeSharpening, terrainWarp } = params;
    const timing = [];

    // Terrain warp — first step, before ocean detection or smoothing
    if (terrainWarp > 0) {
        const t0 = performance.now();
        warpTerrain(mesh, r_elevation, r_xyz, seed, terrainWarp, r_hotspot);
        timing.push({ stage: `Terrain warp (strength=${terrainWarp.toFixed(2)})`, ms: performance.now() - t0 });
    }

    const r_isOcean = new Uint8Array(mesh.numRegions);
    for (let r = 0; r < mesh.numRegions; r++) {
        if (r_elevation[r] <= 0) r_isOcean[r] = 1;
    }

    const preErosion = new Float32Array(r_elevation);

    if (smoothing > 0) {
        const smoothIters = Math.round(1 + smoothing * 4);
        const smoothStr = 0.2 + smoothing * 0.5;
        const t0 = performance.now();
        smoothElevation(mesh, r_elevation, r_isOcean, smoothIters, smoothStr);
        timing.push({ stage: `Smoothing (${smoothIters} iters, str=${smoothStr.toFixed(2)})`, ms: performance.now() - t0 });
    }

    {
        const t0 = performance.now();
        applyDetailNoise(mesh, r_xyz, r_elevation, r_isOcean, seed, {
            dampenField: r_dampen ?? null,
            dampenStrength: DETAIL_NOISE_DAMPEN_STRENGTH,
            amplitudeField: r_orogenic ?? null,
        });
        timing.push({ stage: 'Detail noise L1 (0-100m bumps)', ms: performance.now() - t0 });
    }

    {
        const t0 = performance.now();
        applyDetailNoise(mesh, r_xyz, r_elevation, r_isOcean, seed, {
            amplitudeKm: 0.05,
            frequencyMult: 2.0,
            warpAmpMult: 2.0,
            bipolar: true,
            biasExponent: 0.4,
            seedOffset: 13579,
            dampenField: r_dampen ?? null,
            dampenStrength: DETAIL_NOISE_DAMPEN_STRENGTH,
            amplitudeField: r_orogenic ?? null,
        });
        timing.push({ stage: 'Detail noise L2 (±50m biased)', ms: performance.now() - t0 });
    }

    // ---- Endorheic basins ----
    // Detected on the PRE-CONDITIONING surface: everything shaping-related has
    // run (warp, smoothing, both detail-noise passes) but no drainage
    // conditioning has touched the terrain yet. This is the last moment the
    // depressions are the ones the geology actually produced.
    let basinState = null;
    let basinSelectionSource = null;
    if (basinOpts && basinOpts.enabled !== false) {
        const t0 = performance.now();
        const cellArea = basinOpts.cellArea || null;
        const preConditioningElev = new Float32Array(r_elevation);
        const catalogue = detectBasins(mesh, r_elevation, r_isOcean,
            { cellArea, radiusKm: basinOpts.radiusKm, reliefScale: basinOpts.reliefScale });
        const { selected, warnings, source, carvedByList } = selectBasins(catalogue, basinOpts);
        basinSelectionSource = source;
        attachHypsometry(selected, r_elevation, cellArea, undefined, basinOpts.reliefScale);
        // Cut partial outlets BEFORE building protection, so the floors are
        // computed against the already-notched rim and the invariant holds
        // against what the terrain actually is.
        const incision = inciseOutlets(mesh, r_elevation, selected, { radiusKm: basinOpts.radiusKm });

        const protection = selected.length > 0
            ? buildBasinProtection(mesh, selected, {
                radiusKm: basinOpts.radiusKm,
                ringKm: basinOpts.divideRingKm,
                ringCells: basinOpts.divideRingCells,
              })
            : null;
        basinState = { catalogue, selected, warnings, protection, preConditioningElev, cellArea, source: basinSelectionSource, incision, carvedByList: carvedByList ?? 0 };
        timing.push({
            stage: `Basin detection (${catalogue.length} depressions, ${selected.length} preserved)`,
            ms: performance.now() - t0,
        });
    }
    const protection = basinState ? basinState.protection : null;

    // ---- Lithology ----
    // Classified after basins so closed basins can take evaporites, and before
    // erosion so the landscape can respond to rock strength. Both the rock map
    // and the erodibility field are built from the pre-conditioning surface.
    let lithoState = null;
    if (lithoOpts && lithoOpts.enabled !== false && lithoOpts.tectonics) {
        const t0 = performance.now();
        const litho = classifyLithology(mesh, r_xyz, r_elevation, lithoOpts.tectonics,
            lithoOpts.debugLayers,
            { basins: basinState && basinState.protection
                ? {
                    r_isEndorheic: basinState.protection.isMember,
                    // Zoned salt crust vs clastic fill. Built here because it
                    // needs both the protection's basinIndex and the selected
                    // basins' sink/spill elevations, and computed against the
                    // pre-conditioning surface the basins were catalogued on —
                    // the same surface classifyLithology is reading.
                    r_isSaltCrust: saltCrustMask(
                        basinState.selected, basinState.protection.basinIndex, r_elevation),
                  } : null });
        lithoState = buildLithoState(litho, r_elevation, lithoOpts.cellArea, lithoOpts.strength);
        lithoState.initial = litho;
        timing.push({ stage: 'Lithology classification', ms: performance.now() - t0 });
    }

    if (glacialErosion > 0 || hydraulicErosion > 0 || thermalErosion > 0) {
        const gIters = Math.round(glacialErosion * 10);
        const hIters = Math.round(hydraulicErosion * 20);
        const hK = hydraulicErosion * 0.0006;
        const tIters = Math.round(thermalErosion * 10);
        const talusSlope = 1.2 - thermalErosion * 0.4;
        const kThermal = thermalErosion * 0.15;
        const t0 = performance.now();
        erodeComposite(mesh, r_elevation, r_xyz, r_isOcean,
            hIters, hK, 0.5, 1.0,
            tIters, talusSlope, kThermal,
            gIters, glacialErosion,
            neighborDist,
            protection,
            lithoState,
            hydro);
        timing.push({ stage: `Erosion composite (h=${hIters}, t=${tIters}, g=${gIters})`, ms: performance.now() - t0 });
    }

    if (ridgeSharpening > 0) {
        const rsIters = Math.round(1 + ridgeSharpening * 3);
        const rsStr = ridgeSharpening * 0.08;
        const t0 = performance.now();
        sharpenRidges(mesh, r_elevation, r_isOcean, rsIters, rsStr);
        timing.push({ stage: `Ridge sharpening (${rsIters} iters)`, ms: performance.now() - t0 });
    }

    {
        const t0 = performance.now();
        applySoilCreep(mesh, r_elevation, r_isOcean, 3, 0.1125);
        timing.push({ stage: 'Soil creep (3 iters)', ms: performance.now() - t0 });
    }

    const dl_erosionDelta = new Float32Array(mesh.numRegions);
    for (let r = 0; r < mesh.numRegions; r++) {
        dl_erosionDelta[r] = r_elevation[r] - preErosion[r];
    }

    return { dl_erosionDelta, postTiming: timing, basinState, lithoState };
}

// ─────────────────────────────────────────────────────────────────────────
//  Plate table
//
//  Plate state lives in sparse objects keyed by seed-region ID, which is
//  awkward to serialize and impossible to index from a numeric array. This
//  flattens it into a dense per-plate record list, in a stable order (sorted
//  by seed ID), plus the index remapping needed to interpret r_plate.
// ─────────────────────────────────────────────────────────────────────────
export function buildPlateTable(plateSeeds, plateVec, plateIsOcean, plateDensity, plateDebug,
                                coarseMesh, coarse_xyz, coarse_r_plate) {
    const seedArr = Array.from(plateSeeds).sort((a, b) => a - b);
    const area = {}, centroid = {};
    for (const pid of seedArr) { area[pid] = 0; centroid[pid] = [0, 0, 0]; }
    for (let r = 0; r < coarseMesh.numRegions; r++) {
        const pid = coarse_r_plate[r];
        if (area[pid] === undefined) continue;
        area[pid]++;
        centroid[pid][0] += coarse_xyz[3 * r];
        centroid[pid][1] += coarse_xyz[3 * r + 1];
        centroid[pid][2] += coarse_xyz[3 * r + 2];
    }

    const plates = seedArr.map((pid, index) => {
        const n = area[pid] || 1;
        let [cx, cy, cz] = centroid[pid];
        const clen = Math.sqrt(cx * cx + cy * cy + cz * cz) || 1;
        cx /= clen; cy /= clen; cz /= clen;
        const pv = plateVec[pid] || { pole: [0, 0, 1], omega: 0 };
        const d = plateDebug ? plateDebug[pid] : null;
        return {
            index,
            seedRegion: pid,
            isOcean: plateIsOcean.has(pid),
            density: plateDensity[pid] ?? null,
            // Euler rotation pole (unit vector) and angular velocity. Surface
            // velocity at position p is omega * (pole × p).
            eulerPole: [pv.pole[0], pv.pole[1], pv.pole[2]],
            omega: pv.omega,
            areaFraction: n / coarseMesh.numRegions,
            centroid: [cx, cy, cz],
            centroidLat: Math.asin(Math.max(-1, Math.min(1, cy))) * 180 / Math.PI,
            centroidLon: Math.atan2(cx, cz) * 180 / Math.PI,
            // Plate-physics diagnostics (null when physics didn't run)
            continentalDrag:  d ? d.continentalDrag ?? null : null,
            sizeVelFactor:    d ? d.sizeVelFactor ?? null : null,
            omegaAfter:       d ? d.omegaAfter ?? null : null,
            mantleAlignment:  d ? d.mantleAlignment ?? null : null,
            slabPullStrength: d ? d.slabPullStrength ?? null : null,
            ridgePushStrength: d ? d.ridgePushStrength ?? null : null,
        };
    });

    // r_plate stores seed-region IDs, not 0..numPlates-1. Consumers want dense
    // indices, so ship the lookup rather than making them rebuild it.
    const seedToIndex = new Map(seedArr.map((pid, i) => [pid, i]));
    return { plates, seedArr, seedToIndex };
}

/** Dense 0..numPlates-1 plate index per region, derived from the plate table. */
export function densePlateIndex(r_plate, seedToIndex) {
    const out = new Int32Array(r_plate.length);
    for (let r = 0; r < r_plate.length; r++) {
        const idx = seedToIndex.get(r_plate[r]);
        out[r] = idx === undefined ? -1 : idx;
    }
    return out;
}

// ─────────────────────────────────────────────────────────────────────────
//  Main pipeline
// ─────────────────────────────────────────────────────────────────────────

/**
 * Build a planet from scratch.
 *
 * @param {object} params  slider-equivalent parameters (see defaults below)
 * @param {function} onProgress  (pct, label) => void
 * @returns {object} the full generation context — every intermediate the
 *   pipeline produced, including tectonic internals and the plate table.
 */
export function runGeneratePipeline(params, onProgress = () => {}) {
    const {
        N, P, jitter, nMag, numContinents,
        smoothing, hydraulicErosion, thermalErosion, ridgeSharpening, glacialErosion, terrainWarp,
        continentSizeVariety = 0,
        temperatureOffset = 0, precipitationOffset = 0, landCoverage = 0.3,
        seed: overrideSeed, toggledIndices, skipClimate,
        preserveBasins = true,          // preserve resolvable closed depressions
        preserveBasin = [],             // explicit basin ids, bypass the size floors
        basinMinDepthKm, basinMinAreaKm2, basinMinCells,
        basinDivideRingKm, basinDivideRingCells,
        preserveBasinList = null,       // [{id, retain}] — replaces threshold selection
        carveBasinList = null,          // [{id, retain}] — subtracted from it
        inlandWaterLevels = null,       // { basinId: elevationKm } — supplied, never inferred
        lithology = true,               // classify rock and let erosion respond to it
        lithologyStrength,              // 0 = uniform erodibility, 1 = full rock contrast
        planet: planetSpec = null,      // radius / gravity / rotation; Earth when absent
    } = params;
    const spread = 5;
    const timing = [];
    const tTotal0 = performance.now();

    // Radius sets every physical length scale; gravity scales maximum relief.
    // Earth defaults are exact, so an unnamed planet is bit-identical to a run
    // from before planetary parameters existed.
    const planet = makePlanet(planetSpec || {});

    onProgress(0, 'Shaping the world…');
    const seed = overrideSeed ?? Math.floor(Math.random() * 16777216);
    const rng = makeRng(seed);

    let t0 = performance.now();
    const { mesh, r_xyz } = buildSphere(N, jitter, rng);
    timing.push({ stage: 'Sphere mesh (Fibonacci + Delaunay + pole)', ms: performance.now() - t0 });

    t0 = performance.now();
    const neighborDist = computeNeighborDist(mesh, r_xyz);
    timing.push({ stage: 'Neighbor distances', ms: performance.now() - t0 });

    t0 = performance.now();
    const t_xyz = generateTriangleCenters(mesh, r_xyz);
    timing.push({ stage: 'Triangle centers', ms: performance.now() - t0 });

    onProgress(10, 'Generating coarse plates…');
    t0 = performance.now();
    const { coarseMesh, coarse_xyz, coarse_r_plate, coarsePlateSeeds, coarsePlateVec, coarsePlateIsOcean } =
        generateCoarsePlates(seed, P, numContinents, continentSizeVariety, landCoverage);
    timing.push({ stage: `Coarse plates (${P} plates, ${numContinents} continents)`, ms: performance.now() - t0 });

    onProgress(20, 'Projecting plates…');
    t0 = performance.now();
    const r_plate = projectCoarsePlates(mesh, r_xyz, coarseMesh, coarse_xyz, coarse_r_plate, seed, P);
    timing.push({ stage: 'Project coarse → hi-res', ms: performance.now() - t0 });

    onProgress(25, 'Smoothing boundaries…');
    t0 = performance.now();
    smoothAndReconnectPlates(mesh, r_plate, coarsePlateSeeds, 3);
    timing.push({ stage: 'Smooth projected plates', ms: performance.now() - t0 });

    const plateSeeds = coarsePlateSeeds;
    const plateVec = coarsePlateVec;
    const plateIsOcean = coarsePlateIsOcean;

    // Euler poles before the physics pass, so the exporter can show what the
    // slab-pull / drag / mantle-flow biasing actually changed.
    const plateVecInitial = {};
    for (const pid of plateSeeds) {
        plateVecInitial[pid] = { pole: [...plateVec[pid].pole], omega: plateVec[pid].omega };
    }

    const originalPlateIsOcean = new Set(plateIsOcean);

    if (toggledIndices && toggledIndices.length > 0) {
        const seedArr = Array.from(plateSeeds);
        for (const i of toggledIndices) {
            if (i < seedArr.length) {
                const r = seedArr[i];
                if (plateIsOcean.has(r)) plateIsOcean.delete(r);
                else plateIsOcean.add(r);
            }
        }
    }

    const plateDensity = {};
    const plateDensityLand = {};
    const plateDensityOcean = {};
    for (const r of plateSeeds) {
        const drng = makeRng(r + 777);
        plateDensityOcean[r] = 3.0 + drng() * 0.5;
        plateDensityLand[r] = 2.4 + drng() * 0.5;
        plateDensity[r] = plateIsOcean.has(r) ? plateDensityOcean[r] : plateDensityLand[r];
    }

    const noise = new SimplexNoise(seed);

    // Apply physically-motivated plate motion biasing
    t0 = performance.now();
    const { plateDebug, mantleField, velDelta } = applyPlatePhysics(
        plateVec, plateSeeds, plateIsOcean,
        coarse_r_plate, coarseMesh, coarse_xyz, seed
    );
    timing.push({ stage: 'Plate physics (drag + slab pull + ridge push + mantle flow)', ms: performance.now() - t0 });

    // Build super plates for broad orogenic belts (skip if too few plates)
    let superPlateData = null;
    if (P >= 8) {
        t0 = performance.now();
        superPlateData = buildSuperPlates(coarseMesh, coarse_r_plate, plateSeeds, plateVec, plateIsOcean, plateDensity, r_plate);
        timing.push({ stage: `Super plates (${superPlateData.numSuperPlates} groups from ${P} plates)`, ms: performance.now() - t0 });

        // Apply plate physics to super plates with stronger blending
        t0 = performance.now();
        const spSeeds = new Set();
        for (let i = 0; i < superPlateData.numSuperPlates; i++) spSeeds.add(i);
        applyPlatePhysics(
            superPlateData.superPlateVec, spSeeds, superPlateData.superPlateIsOcean,
            superPlateData.r_superPlate, mesh, r_xyz, seed + 7777,
            SUPER_PLATE_PHYSICS_MULT
        );
        timing.push({ stage: 'Super plate physics', ms: performance.now() - t0 });
    }

    // Expand mantle field from coarse mesh to hi-res via plate averages
    const r_mantleField = new Float32Array(mesh.numRegions);
    {
        const plateMantleSum = {}, plateMantleN = {};
        for (let r = 0; r < coarseMesh.numRegions; r++) {
            const pid = coarse_r_plate[r];
            plateMantleSum[pid] = (plateMantleSum[pid] || 0) + mantleField[r];
            plateMantleN[pid] = (plateMantleN[pid] || 0) + 1;
        }
        for (let r = 0; r < mesh.numRegions; r++) {
            const pid = r_plate[r];
            r_mantleField[r] = plateMantleN[pid] ? plateMantleSum[pid] / plateMantleN[pid] : 0;
        }
    }

    onProgress(35, 'Raising mountains…');
    t0 = performance.now();
    const { r_elevation, mountain_r, coastline_r, ocean_r, r_stress, debugLayers, tectonics, topologyFixup, _timing } =
        assignElevation(mesh, r_xyz, plateIsOcean, r_plate, plateVec, plateSeeds, noise, nMag, seed, spread, plateDensity, superPlateData, r_mantleField,
            { preserveClosedBasins: preserveBasins, radiusKm: planet.radiusKm });
    timing.push({ stage: 'Elevation (collisions + stress + distance fields + assignment)', ms: performance.now() - t0 });

    const prePostElev = new Float32Array(r_elevation);
    const r_dampen = computeDetailDampenField(debugLayers);
    const r_orogenic = computeOrogenicField(debugLayers);

    // True per-cell areas from the dual mesh — basin area, volume and
    // hypsometry are all in km², so a uniform approximation would bias every
    // basin toward whichever latitude its cells happen to sit at.
    const cellArea = regionCellArea(mesh, r_xyz, t_xyz, planet.radiusKm);

    onProgress(60, 'Eroding terrain…');
    t0 = performance.now();
    const { dl_erosionDelta, postTiming, basinState, lithoState } = runPostProcessing(mesh, r_xyz, r_elevation,
        { smoothing, glacialErosion, hydraulicErosion, thermalErosion, ridgeSharpening, terrainWarp },
        neighborDist, seed, debugLayers.hotspot, r_dampen, r_orogenic,
        {
            enabled: preserveBasins,
            cellArea,
            explicitIds: preserveBasin,
            minDepthKm: basinMinDepthKm,
            minAreaKm2: basinMinAreaKm2,
            minCells: basinMinCells,
            radiusKm: planet.radiusKm,
            // Gravity's 1/g relief scaling, so every ...Km the basin stage
            // publishes is on the same vertical scale as `elevation_km`.
            reliefScale: planet.reliefScale,
            divideRingKm: basinDivideRingKm,
            divideRingCells: basinDivideRingCells,
            preserveList: preserveBasinList,
            carveList: carveBasinList,
        },
        {
            enabled: lithology,
            tectonics,
            debugLayers,
            cellArea,
            strength: lithologyStrength,
            radiusKm: planet.radiusKm,
        },
        { cellArea, radiusKm: planet.radiusKm });
    timing.push({ stage: 'Terrain post-processing (total)', ms: performance.now() - t0 });
    debugLayers.erosionDelta = dl_erosionDelta;

    // ---- Basin measurement against the finished terrain ----
    // The catalogue above describes the pre-conditioning surface. Erosion moves
    // divides, so hydrology downstream needs the catchments of the terrain that
    // actually shipped, plus proof the rims survived.
    let basins = null;
    if (basinState) {
        t0 = performance.now();
        const finalIsOcean = new Uint8Array(mesh.numRegions);
        for (let r = 0; r < mesh.numRegions; r++) if (r_elevation[r] <= 0) finalIsOcean[r] = 1;

        const finalMeasurements = remeasureBasins(mesh, r_elevation, finalIsOcean, basinState.selected, cellArea,
            planet.reliefScale);
        // Root the final-terrain catchments at the final-terrain low point.
        // remeasureBasins already found it; using the catalogue sink instead
        // roots the tree on a cell that water flows out of once erosion has
        // moved the basin floor. See computeFinalPreservedCatchments.
        const finalCatchments = computeFinalPreservedCatchments(
            mesh, r_elevation, finalIsOcean, basinState.selected, cellArea,
            { sinks: finalMeasurements.map(f => f.sink) });

        // How far each preserved basin's catchment moved between the
        // pre-conditioning and final terrain.
        const preIsOcean = new Uint8Array(mesh.numRegions);
        for (let r = 0; r < mesh.numRegions; r++) {
            if (basinState.preConditioningElev[r] <= 0) preIsOcean[r] = 1;
        }
        const preCatchments = computeFinalPreservedCatchments(
            mesh, basinState.preConditioningElev, preIsOcean, basinState.selected, cellArea);
        const catchmentDrift = basinState.selected.map((b, i) => ({
            id: b.id,
            ...compareCatchments(
                preCatchments[i].catchmentCells,
                finalCatchments[i].catchmentCells,
                mesh.numRegions),
        }));

        const water = applyInlandWaterLevels(mesh, r_elevation, basinState.selected, inlandWaterLevels || {});
        // Open-ocean connectivity, not elevation sign — a below-sea-level closed
        // basin must not be classed as sea.
        const { isOpenOcean } = findOpenOcean(mesh, finalIsOcean);
        const surfaceClass = surfaceClasses(r_elevation, isOpenOcean, water.isInlandWater);

        // Per-region layers so the basins are reachable from the exporter.
        const r_basinIndex = basinState.protection
            ? basinState.protection.basinIndex
            : new Int32Array(mesh.numRegions).fill(-1);
        const r_isEndorheic = new Uint8Array(mesh.numRegions);
        for (let r = 0; r < mesh.numRegions; r++) if (r_basinIndex[r] >= 0) r_isEndorheic[r] = 1;
        // Distinct from the above: membership is "this cell is basin floor",
        // sink is "this cell is where the water ends up". Conflating them is how
        // flow accumulation came to terminate at basin rims.
        const r_isBasinSink = new Uint8Array(mesh.numRegions);
        for (const f of finalMeasurements) r_isBasinSink[f.sink] = 1;

        // Canonical land for every published fraction: the SUBAERIAL surface.
        // Not `elevation > 0` — that is land_mask, which drops the dry
        // sub-sea-level basin floors and would inflate an endorheic percentage
        // by ~3.5 points on this planet.
        const r_isSubaerial = new Uint8Array(mesh.numRegions);
        let subaerialAreaKm2 = 0;
        for (let r = 0; r < mesh.numRegions; r++) {
            if (surfaceClass[r] === 1) { r_isSubaerial[r] = 1; subaerialAreaKm2 += cellArea[r]; }
        }

        const drainageConsistency = auditCatchmentConsistency(
            finalCatchments[0] ? finalCatchments[0].terminal : new Int32Array(0),
            finalCatchments, cellArea,
            { landAreaKm2: subaerialAreaKm2, landDenominator: 'surface_class == land (subaerial)' });

        basins = {
            catalogue: basinState.catalogue,
            selected: basinState.selected,
            warnings: [...basinState.warnings, ...water.warnings],
            finalMeasurements,
            drainageConsistency,
            sharedCatchmentRoots: finalCatchments.sharedRoots || [],
            finalCatchments: finalCatchments.map(c => ({
                id: c.id, sink: c.sink, catalogueSink: c.catalogueSink,
                catchmentAreaKm2: c.catchmentAreaKm2,
                basinFloorAreaKm2: c.basinFloorAreaKm2,
                catchmentToFloorRatio: c.catchmentToFloorRatio,
            })),
            catchmentDrift,
            waterLevels: water.applied,
            selectionSource: basinState.source,
            incision: basinState.incision,
            carvedByList: basinState.carvedByList,
            drainageHypothesis: {
                preserveListGiven: preserveBasinList ? preserveBasinList.length : 0,
                carveListGiven: carveBasinList ? carveBasinList.length : 0,
                carvedByRetainZero: basinState.carvedByList,
                entries: (preserveBasinList || carveBasinList || []).map(e => ({ id: e.id, retain: e.retain })),
            },
            resolution: basinResolutionContext(mesh.numRegions, planet.radiusKm, {
                minDepthKm: basinMinDepthKm, minAreaKm2: basinMinAreaKm2, minCells: basinMinCells,
            }),
            r_basinIndex,
            r_isEndorheic,
            r_isBasinSink,
            r_isSubaerial,
            subaerialAreaKm2,
            r_surfaceClass: surfaceClass,
            r_inlandWaterBasin: water.inlandWaterBasin,
            r_inlandWaterLevelKm: water.waterLevelKm,
            r_finalTerminal: finalCatchments.length > 0 ? finalCatchments[0].terminal : null,
            preConditioningElev: basinState.preConditioningElev,
            topologyFixup,
        };

        debugLayers.basinIndex = Float32Array.from(r_basinIndex);
        debugLayers.endorheic = Float32Array.from(r_isEndorheic);
        debugLayers.surfaceClass = Float32Array.from(surfaceClass);
        timing.push({ stage: 'Basin remeasurement + final catchments', ms: performance.now() - t0 });
    }

    // Expand plate physics diagnostics to hi-res mesh
    {
        const ppd = expandPlatePhysicsDebug(
            plateDebug, mantleField, velDelta, r_plate, mesh.numRegions,
            coarse_r_plate, coarseMesh.numRegions
        );
        debugLayers.continentalDrag = ppd.dl_continentalDrag;
        debugLayers.sizeVelocity = ppd.dl_sizeVelocity;
        debugLayers.plateSpeed = ppd.dl_plateSpeed;
        debugLayers.velChange = ppd.dl_velChange;
        debugLayers.mantleFlow = ppd.dl_mantleFlow;
    }

    let windResult = null, oceanResult = null, precipResult = null, tempResult = null;

    if (!skipClimate) {
        onProgress(70, 'Simulating wind patterns…');
        t0 = performance.now();
        windResult = computeWind(mesh, r_xyz, r_elevation, plateIsOcean, r_plate, noise);
        timing.push({ stage: 'Wind simulation', ms: performance.now() - t0 });
        if (windResult._windTiming) timing.push(...windResult._windTiming);
        debugLayers.pressureSummer = windResult.r_pressure_summer;
        debugLayers.pressureWinter = windResult.r_pressure_winter;
        debugLayers.windSpeedSummer = windResult.r_wind_speed_summer;
        debugLayers.windSpeedWinter = windResult.r_wind_speed_winter;
        debugLayers.continentality = windResult.r_continentality;

        onProgress(78, 'Computing ocean currents…');
        t0 = performance.now();
        oceanResult = computeOceanCurrents(mesh, r_xyz, r_elevation, windResult);
        timing.push({ stage: 'Ocean currents', ms: performance.now() - t0 });
        if (oceanResult._oceanTiming) timing.push(...oceanResult._oceanTiming);

        onProgress(82, 'Computing precipitation…');
        t0 = performance.now();
        precipResult = computePrecipitation(mesh, r_xyz, r_elevation, windResult, oceanResult, precipitationOffset, landCoverage);
        timing.push({ stage: 'Precipitation', ms: performance.now() - t0 });
        if (precipResult._precipTiming) timing.push(...precipResult._precipTiming);
        debugLayers.precipSummer = precipResult.r_precip_summer;
        debugLayers.precipWinter = precipResult.r_precip_winter;
        debugLayers.rainShadowSummer = precipResult.r_rainshadow_summer;
        debugLayers.rainShadowWinter = precipResult.r_rainshadow_winter;

        onProgress(86, 'Computing temperature…');
        t0 = performance.now();
        tempResult = computeTemperature(mesh, r_xyz, r_elevation, windResult, oceanResult, precipResult, temperatureOffset);
        timing.push({ stage: 'Temperature', ms: performance.now() - t0 });
        if (tempResult._tempTiming) timing.push(...tempResult._tempTiming);
        debugLayers.tempSummer = tempResult.r_temperature_summer;
        debugLayers.tempWinter = tempResult.r_temperature_winter;
        debugLayers.tempContinentality = tempResult.r_tempContinentality;

        t0 = performance.now();
        debugLayers.koppen = classifyKoppen(mesh, r_elevation, tempResult, precipResult);
        timing.push({ stage: 'Köppen classification', ms: performance.now() - t0 });
    }

    onProgress(skipClimate ? 75 : 90, 'Computing triangle elevations…');
    t0 = performance.now();
    const t_elevation = computeTriangleElevations(mesh, r_elevation);
    timing.push({ stage: 'Triangle elevations', ms: performance.now() - t0 });

    // ---- Lithology result, measured against the finished terrain ----
    let lithologyResult = null;
    if (lithoState) {
        const surface = finalSurfaceRock(lithoState);
        const initial = lithoState.initial;
        let exhumed = 0;
        for (let r = 0; r < mesh.numRegions; r++) {
            if (initial.surface[r] !== surface[r]) exhumed++;
        }
        // Where escarpments would form. Computed against the FINAL terrain and
        // the post-erosion cover state, because it is erosion stripping the
        // cover that daylights the interface in the first place.
        const scarp = computeScarpPotential(mesh, r_elevation, lithoState, neighborDist,
            { radiusKm: planet.radiusKm, isLand: basins ? basins.r_isSubaerial : null });
        let scarpCells = 0;
        for (let r = 0; r < mesh.numRegions; r++) if (scarp[r] > 0.25) scarpCells++;

        // Canonical land for every published fraction: the SUBAERIAL surface —
        // not connected to the world ocean and not under an inland lake. Not
        // `elevation > 0`, which is land_mask and silently drops the dry
        // sub-sea-level basin floors. Named in the manifest either way, because
        // a fraction whose denominator is not stated is the trap itself.
        let landDenominator = 'surface_class == land (subaerial)';
        let r_isSubaerial;
        if (basins && basins.r_isSubaerial) {
            r_isSubaerial = basins.r_isSubaerial;
        } else {
            r_isSubaerial = new Uint8Array(mesh.numRegions);
            // No basin pass, so inland water and dry floor are indistinguishable.
            landDenominator = 'elevation > 0 (land_mask) — no basin pass ran';
            for (let r = 0; r < mesh.numRegions; r++) {
                r_isSubaerial[r] = r_elevation[r] > 0 ? 1 : 0;
            }
        }

        lithologyResult = {
            rockClasses: ROCK_CLASSES,
            r_scarpPotential: scarp,
            r_albedo: buildAlbedo(surface),
            scarpCells,
            r_surfaceRock: surface,
            r_basementRock: lithoState.basement,
            r_coverRock: lithoState.cover,
            r_coverThicknessKm: lithoState.coverThicknessKm,
            r_erodibility: lithoState.erodibility,
            r_surfaceRockPreErosion: initial.surface,
            r_coverThicknessPreErosionKm: initial.coverThicknessKm,
            exhumedCells: exhumed,
            erodibilityMeanBeforeNormalisation: lithoState.meanBeforeNormalisation,
            strength: lithologyStrength,
            composition: rockComposition(surface, r_isSubaerial, cellArea),
            compositionSeafloor: rockComposition(surface, null, cellArea),
            landDenominator,
            landAreaKm2: (() => { let a = 0;
                for (let r = 0; r < mesh.numRegions; r++) if (r_isSubaerial[r]) a += cellArea[r];
                return a; })(),
        };
        debugLayers.surfaceRock = Float32Array.from(surface);
        debugLayers.erodibility = Float32Array.from(lithoState.erodibility);
        debugLayers.coverThickness = Float32Array.from(lithoState.coverThicknessKm);
        debugLayers.scarpPotential = Float32Array.from(scarp);
    }

    // ---- Hydrology and ocean topology on the finished terrain ----
    t0 = performance.now();
    const finalIsOceanAll = new Uint8Array(mesh.numRegions);
    for (let r = 0; r < mesh.numRegions; r++) if (r_elevation[r] <= 0) finalIsOceanAll[r] = 1;

    const rivers = computeRivers(mesh, r_elevation, finalIsOceanAll, cellArea, {
        radiusKm: planet.radiusKm,
        // The SINK set, not the member mask: a river entering a preserved basin
        // must run all the way to its low point so the sink's flow accumulation
        // is the inflow a lake water balance needs.
        terminalRoot: basins && basins.r_isBasinSink ? basins.r_isBasinSink : null,
    });
    const oceanBasins = computeOceanBasins(mesh, r_elevation, finalIsOceanAll, cellArea,
        { radiusKm: planet.radiusKm });
    timing.push({ stage: `Rivers + ocean basins (${rivers.mouths.length} mouths, `
        + `${oceanBasins.basins.length} marginal seas)`, ms: performance.now() - t0 });

    const hydrology = {
        r_flowAccumKm2: rivers.flowAccumKm2,
        r_drainTo: rivers.drainTo,
        mouths: rivers.mouths.slice(0, 200),
        mouthCount: rivers.mouths.length,
        dischargeDistribution: [1e3, 1e4, 1e5, 1e6].map(t => ({
            aboveKm2: t,
            outlets: rivers.mouths.reduce((n, m) => n + (m.dischargeAreaKm2 > t ? 1 : 0), 0),
        })),
        oceanBasins: oceanBasins.basins.map(b => {
            const { members, ...rest } = b;
            return rest;
        }),
        r_oceanBasinIndex: oceanBasins.label,
    };
    debugLayers.flowAccum = Float32Array.from(rivers.flowAccumKm2);
    debugLayers.oceanBasin = Float32Array.from(oceanBasins.label);

    const plateTable = buildPlateTable(plateSeeds, plateVec, plateIsOcean, plateDensity, plateDebug,
        coarseMesh, coarse_xyz, coarse_r_plate);

    return {
        mesh, r_xyz, t_xyz, neighborDist,
        r_plate, plateSeeds, plateVec, plateVecInitial, plateIsOcean, originalPlateIsOcean,
        plateDensity, plateDensityLand, plateDensityOcean,
        plateDebug, plateTable,
        coarseMesh, coarse_xyz, coarse_r_plate,
        superPlateData, mantleField, velDelta, r_mantleField,
        r_elevation, t_elevation, prePostElev,
        mountain_r, coastline_r, ocean_r, r_stress,
        debugLayers, tectonics, basins, lithology: lithologyResult, hydrology, cellArea,
        planet, planetSummary: planetSummary(planet),
        r_dampen, r_orogenic,
        windResult, oceanResult, precipResult, tempResult,
        noise, seed, nMag, P,
        // Every parameter that changes the terrain must appear here, because
        // this is what the export manifest records and therefore the only
        // account a build gives of itself. `lithologyStrength` and `lithology`
        // were missing: a build made with --lithology-strength 0.682 got a
        // different finalElevation hash, so the terrain was distinguishable,
        // but nothing anywhere said WHY it differed. A hash that changes for an
        // unrecorded reason is worse than no hash.
        params: { N, P, jitter, nMag, numContinents, smoothing, terrainWarp, hydraulicErosion,
                  thermalErosion, ridgeSharpening, glacialErosion, continentSizeVariety,
                  temperatureOffset, precipitationOffset, landCoverage, seed,
                  lithology: !!lithology,
                  lithologyStrength: lithology ? (lithologyStrength ?? 1) : null },
        skipClimate: !!skipClimate,
        _timing,
        _pipelineTiming: timing,
        _postTiming: postTiming,
        _workerTotal: performance.now() - tTotal0,
    };
}
