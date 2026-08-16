// Web Worker — runs the pure computation pipeline off the main thread.
// Handles: generate, reapply, editRecompute commands.

import { makeRng } from './rng.js';
import { SimplexNoise } from './simplex-noise.js';
import { setDelaunator, buildSphere, generateTriangleCenters, computeNeighborDist } from './sphere-mesh.js';
import { assignElevation } from './elevation.js';
import { buildSuperPlates } from './super-plates.js';
import { computeWind } from './wind.js';
import { computeOceanCurrents } from './ocean.js';
import { computePrecipitation } from './precipitation.js';
import { computeTemperature } from './temperature.js';
import { classifyKoppen } from './koppen.js';
import { computeTerrainMetrics } from './terrain-metrics.js';
import {
    runGeneratePipeline, computeTriangleElevations, computeDetailDampenField,
    computeOrogenicField, runPostProcessing, densePlateIndex,
} from './pipeline.js';
import Delaunator from 'https://cdn.jsdelivr.net/npm/delaunator@5.0.1/+esm';

setDelaunator(Delaunator);

// Retained state between commands (avoids re-sending mesh for reapply/edit)
let W = null;

function progress(pct, label) {
    self.postMessage({ type: 'progress', pct, label });
}

// Basin options for the paths that re-run erosion (reapply, editRecompute).
// They must rebuild protection from scratch: erosion is what breaches rims, so
// re-running it without protection would undo preservation.
// Lithology options for the paths that re-run erosion. Like basins, this has to
// be rebuilt from scratch: erosion is what responds to rock strength, so
// re-running it without the erodibility field would silently drop the effect.
function lithoOptsFrom(data, tectonics, debugLayers) {
    if (data && data.lithology === false) return null;
    return {
        enabled: true,
        tectonics,
        debugLayers,
        cellArea: W && W.cellArea ? W.cellArea : null,
        strength: data ? data.lithologyStrength : undefined,
    };
}

function basinOptsFrom(data) {
    if (data && data.preserveBasins === false) return null;
    return {
        enabled: true,
        cellArea: W && W.cellArea ? W.cellArea : null,
        explicitIds: (data && data.preserveBasin) || [],
        minDepthKm: data ? data.basinMinDepthKm : undefined,
        minAreaKm2: data ? data.basinMinAreaKm2 : undefined,
        minCells: data ? data.basinMinCells : undefined,
    };
}

function getClimateParams(data) {
    const temperatureOffset = data?.temperatureOffset ?? W?.temperatureOffset ?? 0;
    const precipitationOffset = data?.precipitationOffset ?? W?.precipitationOffset ?? 0;
    const landCoverage = data?.landCoverage ?? W?.landCoverage ?? 0.3;
    if (W) { W.temperatureOffset = temperatureOffset; W.precipitationOffset = precipitationOffset; W.landCoverage = landCoverage; }
    return { temperatureOffset, precipitationOffset, landCoverage };
}

function buildClimateFields(windResult, oceanResult, precipResult, tempResult) {
    return {
        r_wind_east_summer: windResult?.r_wind_east_summer ?? null,
        r_wind_north_summer: windResult?.r_wind_north_summer ?? null,
        r_wind_east_winter: windResult?.r_wind_east_winter ?? null,
        r_wind_north_winter: windResult?.r_wind_north_winter ?? null,
        itczLons: windResult?.itczLons ?? null,
        itczLatsSummer: windResult?.itczLatsSummer ?? null,
        itczLatsWinter: windResult?.itczLatsWinter ?? null,
        r_ocean_current_east_summer: oceanResult?.r_ocean_current_east_summer ?? null,
        r_ocean_current_north_summer: oceanResult?.r_ocean_current_north_summer ?? null,
        r_ocean_current_east_winter: oceanResult?.r_ocean_current_east_winter ?? null,
        r_ocean_current_north_winter: oceanResult?.r_ocean_current_north_winter ?? null,
        r_ocean_speed_summer: oceanResult?.r_ocean_speed_summer ?? null,
        r_ocean_speed_winter: oceanResult?.r_ocean_speed_winter ?? null,
        r_ocean_warmth_summer: oceanResult?.r_ocean_warmth_summer ?? null,
        r_ocean_warmth_winter: oceanResult?.r_ocean_warmth_winter ?? null,
        r_precip_summer: precipResult?.r_precip_summer ?? null,
        r_precip_winter: precipResult?.r_precip_winter ?? null,
        r_temperature_summer: tempResult?.r_temperature_summer ?? null,
        r_temperature_winter: tempResult?.r_temperature_winter ?? null,
    };
}

/**
 * Basin state for the main thread. Member region lists are dropped — they can
 * be millions of indices and the per-region basin_index layer carries the same
 * information in a form the renderer and exporter can both use directly.
 */
function serialiseBasins(b) {
    const light = (x) => { const { members, hypsometry, ...rest } = x; return { ...rest, hypsometry }; };
    return {
        catalogue: b.catalogue.map(light),
        selected: b.selected.map(light),
        warnings: b.warnings,
        finalMeasurements: b.finalMeasurements,
        finalCatchments: b.finalCatchments,
        catchmentDrift: b.catchmentDrift,
        waterLevels: b.waterLevels,
        resolution: b.resolution,
        topologyFixup: b.topologyFixup,
        r_basinIndex: b.r_basinIndex,
        r_isEndorheic: b.r_isEndorheic,
        r_surfaceClass: b.r_surfaceClass,
        r_inlandWaterBasin: b.r_inlandWaterBasin,
        r_inlandWaterLevelKm: b.r_inlandWaterLevelKm,
        preConditioningElev: b.preConditioningElev,
    };
}

function handleGenerate(data) {
    try {
        const ctx = runGeneratePipeline(data, progress);

        const {
            mesh, r_xyz, t_xyz, neighborDist,
            r_plate, plateSeeds, plateVec, plateIsOcean, originalPlateIsOcean,
            plateDensity, plateDensityLand, plateDensityOcean,
            coarseMesh, coarse_r_plate, superPlateData,
            r_elevation, t_elevation, prePostElev,
            mountain_r, coastline_r, ocean_r, r_stress,
            debugLayers, tectonics, plateTable,
            r_dampen, r_orogenic,
            windResult, oceanResult, precipResult, tempResult,
            noise, seed, skipClimate,
        } = ctx;
        const { nMag, P, temperatureOffset, precipitationOffset, landCoverage } = ctx.params;

        let t0 = performance.now();
        // Retain state for reapply/edit (clone what we'll transfer)
        W = {
            mesh, r_xyz: new Float32Array(r_xyz), t_xyz: new Float32Array(t_xyz),
            neighborDist,
            r_plate: new Int32Array(r_plate), plateSeeds: new Set(plateSeeds), plateVec,
            plateIsOcean: new Set(plateIsOcean), originalPlateIsOcean: new Set(originalPlateIsOcean),
            plateDensity: Object.assign({}, plateDensity),
            plateDensityLand: Object.assign({}, plateDensityLand),
            plateDensityOcean: Object.assign({}, plateDensityOcean),
            prePostElev: new Float32Array(prePostElev),
            r_elevation_final: new Float32Array(r_elevation),
            seed, nMag, noise, P,
            mountain_r: new Set(mountain_r), coastline_r: new Set(coastline_r), ocean_r: new Set(ocean_r),
            r_stress: new Float32Array(r_stress),
            temperatureOffset, precipitationOffset, landCoverage,
            cachedWind: windResult, cachedOcean: oceanResult,
            // Retain detail-noise dampen + orogenic fields so reapply (which
            // reuses prePostElev) shapes the noise the same way as the initial
            // generate over craton/basin and orogenic regions.
            r_dampen: r_dampen ? new Float32Array(r_dampen) : null,
            r_orogenic: r_orogenic ? new Float32Array(r_orogenic) : null,
            // Retain coarse-plate data so editRecompute can rebuild super
            // plates with the same detail-independent adjacency graph.
            coarseMesh, coarse_r_plate: new Int32Array(coarse_r_plate),
            cellArea: ctx.cellArea,
        };
        ctx._pipelineTiming.push({ stage: 'Clone state for retention', ms: performance.now() - t0 });

        // Compute terrain quality metrics using retained-state clones
        // (the originals will be transferred and neutered below).
        let terrainMetrics = null;
        try {
            terrainMetrics = computeTerrainMetrics({
                mesh: W.mesh,
                r_xyz: W.r_xyz,
                r_elevation: W.r_elevation_final,
                r_plate: W.r_plate,
                plateIsOcean: Array.from(W.plateIsOcean),
                r_stress: W.r_stress,
                debugLayers,
                prePostElev: W.prePostElev,
            });
        } catch (e) {
            terrainMetrics = { _error: e.message };
        }

        // Build result — typed arrays we no longer need are transferred (zero-copy).
        // mesh.triangles/halfedges are NOT transferred because W.mesh retains them.
        const result = {
            type: 'done',
            triangles: mesh.triangles,
            halfedges: mesh.halfedges,
            numRegions: mesh.numRegions,
            r_xyz, t_xyz, r_plate,
            plateSeeds: Array.from(plateSeeds),
            plateVec,
            plateIsOcean: Array.from(plateIsOcean),
            originalPlateIsOcean: Array.from(originalPlateIsOcean),
            plateDensity, plateDensityLand, plateDensityOcean,
            prePostElev,
            r_elevation, t_elevation,
            mountain_r: Array.from(mountain_r),
            coastline_r: Array.from(coastline_r),
            ocean_r: Array.from(ocean_r),
            r_stress,
            ...buildClimateFields(windResult, oceanResult, precipResult, tempResult),
            skipClimate: !!skipClimate,
            seed, nMag,
            debugLayers,
            // Full tectonic internals — see buildTectonicsBundle in elevation.js.
            tectonics,
            basins: ctx.basins ? serialiseBasins(ctx.basins) : null,
            lithology: ctx.lithology || null,
            plates: plateTable.plates,
            r_plateIndex: densePlateIndex(r_plate, plateTable.seedToIndex),
            superPlates: superPlateData ? {
                r_superPlate: superPlateData.r_superPlate,
                numSuperPlates: superPlateData.numSuperPlates,
                superPlateIsOcean: Array.from(superPlateData.superPlateIsOcean),
                superPlateVec: superPlateData.superPlateVec,
                superPlateDensity: superPlateData.superPlateDensity,
            } : null,
            _timing: ctx._timing,             // elevation sub-stages from assignElevation
            _pipelineTiming: ctx._pipelineTiming,
            _postTiming: ctx._postTiming,
            _workerTotal: ctx._workerTotal,
            _params: ctx.params,
            terrainMetrics
        };

        // Transfer arrays the worker no longer needs (cloned copies kept in W).
        // The tectonics bundle is ~40 full-resolution arrays; at high detail
        // structured-cloning it instead would double peak memory. Nothing in W
        // aliases them — W holds its own copies — so transferring is safe.
        // A buffer must appear at most once, and tectonics.r_stress is the same
        // array as r_stress, hence the Set.
        const buffers = new Set([
            r_xyz.buffer, t_xyz.buffer, r_plate.buffer,
            prePostElev.buffer, r_elevation.buffer, t_elevation.buffer,
            r_stress.buffer
        ]);
        for (const v of Object.values(tectonics)) {
            if (ArrayBuffer.isView(v)) buffers.add(v.buffer);
        }
        if (result.r_plateIndex) buffers.add(result.r_plateIndex.buffer);

        self.postMessage(result, Array.from(buffers));

    } catch (err) {
        self.postMessage({ type: 'error', message: err.message, stack: err.stack });
    }
}
function handleReapply(data) {
    if (!W) { self.postMessage({ type: 'error', message: 'No retained state for reapply' }); return; }

    const skipClimate = !!data.skipClimate;
    const { temperatureOffset, precipitationOffset, landCoverage } = getClimateParams(data);

    try {
        const tTotal0 = performance.now();

        progress(0, 'Reapplying terrain\u2026');

        let t0 = performance.now();
        const r_elevation = new Float32Array(W.prePostElev);
        const tClone = performance.now() - t0;

        progress(20, 'Eroding terrain\u2026');
        t0 = performance.now();
        const { dl_erosionDelta, postTiming } = runPostProcessing(W.mesh, W.r_xyz, r_elevation, data, W.neighborDist, W.seed, undefined, W.r_dampen, W.r_orogenic, basinOptsFrom(data));
        const tPost = performance.now() - t0;

        // Update retained final elevation for deferred climate
        W.r_elevation_final = new Float32Array(r_elevation);

        let windResult = null, oceanResult = null, precipResult = null, tempResult = null;
        let tWind = 0, tOcean = 0, tPrecip = 0, tTemp = 0;

        if (!skipClimate) {
            progress(60, 'Simulating wind patterns\u2026');
            t0 = performance.now();
            windResult = computeWind(W.mesh, W.r_xyz, r_elevation, W.plateIsOcean, W.r_plate, W.noise);
            tWind = performance.now() - t0;

            progress(75, 'Computing ocean currents\u2026');
            t0 = performance.now();
            oceanResult = computeOceanCurrents(W.mesh, W.r_xyz, r_elevation, windResult);
            tOcean = performance.now() - t0;

            progress(80, 'Computing precipitation\u2026');
            t0 = performance.now();
            precipResult = computePrecipitation(W.mesh, W.r_xyz, r_elevation, windResult, oceanResult, precipitationOffset, landCoverage);
            tPrecip = performance.now() - t0;

            progress(85, 'Computing temperature\u2026');
            t0 = performance.now();
            tempResult = computeTemperature(W.mesh, W.r_xyz, r_elevation, windResult, oceanResult, precipResult, temperatureOffset);
            tTemp = performance.now() - t0;

            W.cachedWind = windResult;
            W.cachedOcean = oceanResult;
        } else {
            W.cachedWind = null;
            W.cachedOcean = null;
        }

        progress(skipClimate ? 70 : 90, 'Computing triangle elevations\u2026');
        t0 = performance.now();
        const t_elevation = computeTriangleElevations(W.mesh, r_elevation);
        const tTriElev = performance.now() - t0;

        const tWorkerTotal = performance.now() - tTotal0;

        const result = {
            type: 'reapplyDone',
            skipClimate,
            r_elevation,
            t_elevation,
            erosionDelta: dl_erosionDelta,
            ...buildClimateFields(windResult, oceanResult, precipResult, tempResult),
            windDebugLayers: windResult ? {
                pressureSummer: windResult.r_pressure_summer,
                pressureWinter: windResult.r_pressure_winter,
                windSpeedSummer: windResult.r_wind_speed_summer,
                windSpeedWinter: windResult.r_wind_speed_winter,
                precipSummer: precipResult.r_precip_summer,
                precipWinter: precipResult.r_precip_winter,
                rainShadowSummer: precipResult.r_rainshadow_summer,
                rainShadowWinter: precipResult.r_rainshadow_winter,
                tempSummer: tempResult.r_temperature_summer,
                tempWinter: tempResult.r_temperature_winter,
                koppen: classifyKoppen(W.mesh, r_elevation, tempResult, precipResult)
            } : null,
            _reapplyTiming: {
                clone: tClone,
                postProcessing: tPost,
                wind: tWind,
                ocean: tOcean,
                precipitation: tPrecip,
                temperature: tTemp,
                triangleElevations: tTriElev,
                workerTotal: tWorkerTotal
            },
            _postTiming: postTiming
        };

        self.postMessage(result, [r_elevation.buffer, t_elevation.buffer, dl_erosionDelta.buffer]);

    } catch (err) {
        self.postMessage({ type: 'error', message: err.message, stack: err.stack });
    }
}

function handleEditRecompute(data) {
    if (!W) { self.postMessage({ type: 'error', message: 'No retained state for editRecompute' }); return; }

    const skipClimate = !!data.skipClimate;
    const { temperatureOffset, precipitationOffset, landCoverage } = getClimateParams(data);

    try {
        const tTotal0 = performance.now();

        progress(0, 'Rebuilding elevation\u2026');

        // Update retained plate state
        W.plateIsOcean = new Set(data.plateIsOcean);
        W.plateDensity = Object.assign({}, data.plateDensity);

        const { mesh, r_xyz, plateIsOcean, r_plate, plateVec, plateSeeds, noise, seed } = W;
        const nMag = data.nMag;
        const spread = 5;

        // Rebuild super plates from updated plate ocean/density state
        // (uses retained coarse-plate data for detail-stable adjacency graph)
        let superPlateData = null;
        if ((W.P || 0) >= 8) {
            superPlateData = buildSuperPlates(W.coarseMesh, W.coarse_r_plate, plateSeeds, plateVec, plateIsOcean, W.plateDensity, r_plate);
        }

        let t0 = performance.now();
        const { r_elevation, mountain_r, coastline_r, ocean_r, r_stress, debugLayers, tectonics, _timing } =
            assignElevation(mesh, r_xyz, plateIsOcean, r_plate, plateVec, plateSeeds, noise, nMag, seed, spread, W.plateDensity, superPlateData,
                { preserveClosedBasins: data.preserveBasins !== false });
        const tElev = performance.now() - t0;

        const prePostElev = new Float32Array(r_elevation);
        const r_dampen = computeDetailDampenField(debugLayers);
        const r_orogenic = computeOrogenicField(debugLayers);
        W.r_dampen = r_dampen ? new Float32Array(r_dampen) : null;
        W.r_orogenic = r_orogenic ? new Float32Array(r_orogenic) : null;

        progress(50, 'Eroding terrain\u2026');
        t0 = performance.now();
        const { dl_erosionDelta, postTiming } = runPostProcessing(mesh, r_xyz, r_elevation, data, W.neighborDist, W.seed, debugLayers.hotspot, r_dampen, r_orogenic, basinOptsFrom(data), lithoOptsFrom(data, tectonics, debugLayers));
        const tPost = performance.now() - t0;
        debugLayers.erosionDelta = dl_erosionDelta;

        // Update retained final elevation for deferred climate
        W.r_elevation_final = new Float32Array(r_elevation);

        let windResult = null, oceanResult = null, precipResult = null, tempResult = null;
        let tWind = 0, tOcean = 0, tPrecip = 0, tTemp = 0;

        if (!skipClimate) {
            progress(65, 'Simulating wind patterns\u2026');
            t0 = performance.now();
            windResult = computeWind(mesh, r_xyz, r_elevation, plateIsOcean, r_plate, W.noise);
            tWind = performance.now() - t0;
            debugLayers.pressureSummer = windResult.r_pressure_summer;
            debugLayers.pressureWinter = windResult.r_pressure_winter;
            debugLayers.windSpeedSummer = windResult.r_wind_speed_summer;
            debugLayers.windSpeedWinter = windResult.r_wind_speed_winter;
            debugLayers.continentality = windResult.r_continentality;

            progress(78, 'Computing ocean currents\u2026');
            t0 = performance.now();
            oceanResult = computeOceanCurrents(mesh, r_xyz, r_elevation, windResult);
            tOcean = performance.now() - t0;

            progress(82, 'Computing precipitation\u2026');
            t0 = performance.now();
            precipResult = computePrecipitation(mesh, r_xyz, r_elevation, windResult, oceanResult, precipitationOffset, landCoverage);
            tPrecip = performance.now() - t0;
            debugLayers.precipSummer = precipResult.r_precip_summer;
            debugLayers.precipWinter = precipResult.r_precip_winter;
            debugLayers.rainShadowSummer = precipResult.r_rainshadow_summer;
            debugLayers.rainShadowWinter = precipResult.r_rainshadow_winter;

            progress(86, 'Computing temperature\u2026');
            t0 = performance.now();
            tempResult = computeTemperature(mesh, r_xyz, r_elevation, windResult, oceanResult, precipResult, temperatureOffset);
            tTemp = performance.now() - t0;
            debugLayers.tempSummer = tempResult.r_temperature_summer;
            debugLayers.tempWinter = tempResult.r_temperature_winter;
            debugLayers.tempContinentality = tempResult.r_tempContinentality;

            debugLayers.koppen = classifyKoppen(mesh, r_elevation, tempResult, precipResult);

            W.cachedWind = windResult;
            W.cachedOcean = oceanResult;
        } else {
            W.cachedWind = null;
            W.cachedOcean = null;
        }

        progress(skipClimate ? 75 : 90, 'Computing triangle elevations\u2026');
        t0 = performance.now();
        const t_elevation = computeTriangleElevations(mesh, r_elevation);
        const tTriElev = performance.now() - t0;

        // Update retained state
        t0 = performance.now();
        W.prePostElev = new Float32Array(prePostElev);
        W.mountain_r = new Set(mountain_r);
        W.coastline_r = new Set(coastline_r);
        W.ocean_r = new Set(ocean_r);
        W.r_stress = new Float32Array(r_stress);
        const tRetain = performance.now() - t0;

        const tWorkerTotal = performance.now() - tTotal0;

        const result = {
            type: 'editDone',
            skipClimate,
            prePostElev,
            r_elevation,
            t_elevation,
            mountain_r: Array.from(mountain_r),
            coastline_r: Array.from(coastline_r),
            ocean_r: Array.from(ocean_r),
            r_stress,
            ...buildClimateFields(windResult, oceanResult, precipResult, tempResult),
            debugLayers,
            tectonics,
            superPlates: superPlateData ? {
                r_superPlate: superPlateData.r_superPlate,
                numSuperPlates: superPlateData.numSuperPlates,
                superPlateIsOcean: Array.from(superPlateData.superPlateIsOcean),
                superPlateVec: superPlateData.superPlateVec,
                superPlateDensity: superPlateData.superPlateDensity,
            } : null,
            _editTiming: {
                elevation: tElev,
                postProcessing: tPost,
                wind: tWind,
                ocean: tOcean,
                precipitation: tPrecip,
                temperature: tTemp,
                triangleElevations: tTriElev,
                retainState: tRetain,
                workerTotal: tWorkerTotal
            },
            _timing,        // elevation sub-stages
            _postTiming: postTiming
        };

        self.postMessage(result, [
            prePostElev.buffer, r_elevation.buffer, t_elevation.buffer, r_stress.buffer
        ]);

    } catch (err) {
        self.postMessage({ type: 'error', message: err.message, stack: err.stack });
    }
}

function handleComputeClimate(data) {
    if (!W) { self.postMessage({ type: 'error', message: 'No retained state for computeClimate' }); return; }

    const { temperatureOffset, precipitationOffset, landCoverage } = getClimateParams(data);

    try {
        const tTotal0 = performance.now();
        const { mesh, r_xyz, r_elevation_final, plateIsOcean, r_plate, noise } = W;

        let windResult = W.cachedWind;
        let oceanResult = W.cachedOcean;
        let tWind = 0, tOcean = 0;
        let t0;

        if (!windResult) {
            progress(0, 'Simulating wind patterns\u2026');
            t0 = performance.now();
            windResult = computeWind(mesh, r_xyz, r_elevation_final, plateIsOcean, r_plate, noise);
            tWind = performance.now() - t0;

            progress(30, 'Computing ocean currents\u2026');
            t0 = performance.now();
            oceanResult = computeOceanCurrents(mesh, r_xyz, r_elevation_final, windResult);
            tOcean = performance.now() - t0;

            W.cachedWind = windResult;
            W.cachedOcean = oceanResult;
        }

        progress(50, 'Computing precipitation\u2026');
        t0 = performance.now();
        const precipResult = computePrecipitation(mesh, r_xyz, r_elevation_final, windResult, oceanResult, precipitationOffset, landCoverage);
        const tPrecip = performance.now() - t0;

        progress(70, 'Computing temperature\u2026');
        t0 = performance.now();
        const tempResult = computeTemperature(mesh, r_xyz, r_elevation_final, windResult, oceanResult, precipResult, temperatureOffset);
        const tTemp = performance.now() - t0;

        progress(88, 'Classifying climates\u2026');
        t0 = performance.now();
        const koppen = classifyKoppen(mesh, r_elevation_final, tempResult, precipResult);
        const tKoppen = performance.now() - t0;

        const tWorkerTotal = performance.now() - tTotal0;

        const climateDebugLayers = {
            pressureSummer: windResult.r_pressure_summer,
            pressureWinter: windResult.r_pressure_winter,
            windSpeedSummer: windResult.r_wind_speed_summer,
            windSpeedWinter: windResult.r_wind_speed_winter,
            continentality: windResult.r_continentality,
            precipSummer: precipResult.r_precip_summer,
            precipWinter: precipResult.r_precip_winter,
            rainShadowSummer: precipResult.r_rainshadow_summer,
            rainShadowWinter: precipResult.r_rainshadow_winter,
            tempSummer: tempResult.r_temperature_summer,
            tempWinter: tempResult.r_temperature_winter,
            koppen
        };

        progress(95, 'Done');

        self.postMessage({
            type: 'climateDone',
            r_wind_east_summer: windResult.r_wind_east_summer,
            r_wind_north_summer: windResult.r_wind_north_summer,
            r_wind_east_winter: windResult.r_wind_east_winter,
            r_wind_north_winter: windResult.r_wind_north_winter,
            itczLons: windResult.itczLons,
            itczLatsSummer: windResult.itczLatsSummer,
            itczLatsWinter: windResult.itczLatsWinter,
            r_ocean_current_east_summer: oceanResult.r_ocean_current_east_summer,
            r_ocean_current_north_summer: oceanResult.r_ocean_current_north_summer,
            r_ocean_current_east_winter: oceanResult.r_ocean_current_east_winter,
            r_ocean_current_north_winter: oceanResult.r_ocean_current_north_winter,
            r_ocean_speed_summer: oceanResult.r_ocean_speed_summer,
            r_ocean_speed_winter: oceanResult.r_ocean_speed_winter,
            r_ocean_warmth_summer: oceanResult.r_ocean_warmth_summer,
            r_ocean_warmth_winter: oceanResult.r_ocean_warmth_winter,
            r_precip_summer: precipResult.r_precip_summer,
            r_precip_winter: precipResult.r_precip_winter,
            r_temperature_summer: tempResult.r_temperature_summer,
            r_temperature_winter: tempResult.r_temperature_winter,
            climateDebugLayers,
            _climateTiming: {
                wind: tWind,
                ocean: tOcean,
                precipitation: tPrecip,
                temperature: tTemp,
                koppen: tKoppen,
                workerTotal: tWorkerTotal
            }
        });

    } catch (err) {
        self.postMessage({ type: 'error', message: err.message, stack: err.stack });
    }
}

// ─── Heightmap import ───────────────────────────────────────────────

/** Bilinear interpolation with equirectangular wrapping. */
function sampleBilinear(pixels, imgW, imgH, px, py) {
    // Clamp vertically, wrap horizontally
    py = Math.max(0, Math.min(py, imgH - 1));
    const x0 = Math.floor(px), y0 = Math.floor(py);
    const x1 = (x0 + 1) % imgW;     // horizontal wrap
    const y1 = Math.min(y0 + 1, imgH - 1); // vertical clamp
    const fx = px - x0, fy = py - y0;
    const v00 = pixels[y0 * imgW + ((x0 % imgW) + imgW) % imgW];
    const v10 = pixels[y0 * imgW + x1];
    const v01 = pixels[y1 * imgW + ((x0 % imgW) + imgW) % imgW];
    const v11 = pixels[y1 * imgW + x1];
    return (v00 * (1 - fx) * (1 - fy) +
            v10 * fx * (1 - fy) +
            v01 * (1 - fx) * fy +
            v11 * fx * fy);
}

/**
 * Convert grayscale 0–255 to internal elevation.
 * 0 → -0.5 (ocean floor)
 * 1–255 → inverse of 6·t² so grayscale maps linearly to km.
 * Simple sqrt inversion: t = sqrt((v-1) / 254).
 */
function grayscaleToElevation(v) {
    if (v < 1) return -0.5; // ocean (black pixels; catches interpolated fractional values too)
    return Math.sqrt((v - 1) / 254);
}

/**
 * Sample an equirectangular grayscale heightmap onto sphere mesh regions.
 * Returns r_elevation (Float32Array).
 */
function sampleHeightmap(mesh, r_xyz, imageData, imgW, imgH) {
    const r_elevation = new Float32Array(mesh.numRegions);
    for (let r = 0; r < mesh.numRegions; r++) {
        const x = r_xyz[3 * r], y = r_xyz[3 * r + 1], z = r_xyz[3 * r + 2];
        const lat = Math.asin(Math.max(-1, Math.min(1, y)));
        const lon = Math.atan2(x, z);
        // Map lat/lon → pixel coords (equirectangular)
        const px = (lon / Math.PI + 1) * 0.5 * imgW; // 0..W
        const py = (0.5 - lat / Math.PI) * imgH;     // 0..H
        const gray = sampleBilinear(imageData, imgW, imgH, px, py);
        r_elevation[r] = grayscaleToElevation(gray);
    }
    return r_elevation;
}

/**
 * BFS flood fill to derive synthetic plates from elevation.
 * Creates one "plate" per connected land mass and one per connected ocean basin.
 */
function deriveSyntheticPlates(mesh, r_elevation) {
    const N = mesh.numRegions;
    const r_plate = new Int32Array(N).fill(-1);
    const plateSeeds = new Set();
    const plateIsOcean = new Set();
    const plateVec = {};
    const { adjOffset, adjList } = mesh;

    let plateId = 0;
    for (let r = 0; r < N; r++) {
        if (r_plate[r] >= 0) continue;
        const isOcean = r_elevation[r] <= 0;
        // BFS from this region
        r_plate[r] = r; // use r as the plate seed
        plateSeeds.add(r);
        plateVec[r] = [0, 0, 0]; // zero velocity
        if (isOcean) plateIsOcean.add(r);
        const queue = [r];
        let head = 0;
        while (head < queue.length) {
            const cur = queue[head++];
            const end = adjOffset[cur + 1];
            for (let ni = adjOffset[cur]; ni < end; ni++) {
                const nb = adjList[ni];
                if (r_plate[nb] >= 0) continue;
                const nbOcean = r_elevation[nb] <= 0;
                if (nbOcean === isOcean) {
                    r_plate[nb] = r;
                    queue.push(nb);
                }
            }
        }
        plateId++;
    }

    return { r_plate, plateSeeds, plateIsOcean, plateVec };
}

function handleImportHeightmap(data) {
    const { N, jitter, grayscale, imageWidth, imageHeight, smoothing, hydraulicErosion, thermalErosion, ridgeSharpening, glacialErosion, terrainWarp, temperatureOffset = 0, precipitationOffset = 0, landCoverage = 0.3, seed: overrideSeed, skipClimate } = data;
    const timing = [];

    try {
        const tTotal0 = performance.now();

        progress(0, 'Building sphere mesh\u2026');
        const seed = overrideSeed ?? Math.floor(Math.random() * 16777216);
        const rng = makeRng(seed);

        let t0 = performance.now();
        const { mesh, r_xyz } = buildSphere(N, jitter, rng);
        timing.push({ stage: 'Sphere mesh', ms: performance.now() - t0 });

        t0 = performance.now();
        const neighborDist = computeNeighborDist(mesh, r_xyz);
        timing.push({ stage: 'Neighbor distances', ms: performance.now() - t0 });

        t0 = performance.now();
        const t_xyz = generateTriangleCenters(mesh, r_xyz);
        timing.push({ stage: 'Triangle centers', ms: performance.now() - t0 });

        progress(20, 'Sampling heightmap\u2026');
        t0 = performance.now();
        const r_elevation = sampleHeightmap(mesh, r_xyz, grayscale, imageWidth, imageHeight);
        timing.push({ stage: 'Sample heightmap', ms: performance.now() - t0 });

        const prePostElev = new Float32Array(r_elevation);

        progress(35, 'Processing terrain\u2026');
        t0 = performance.now();
        const { dl_erosionDelta, postTiming } = runPostProcessing(mesh, r_xyz, r_elevation, { smoothing, glacialErosion, hydraulicErosion, thermalErosion, ridgeSharpening, terrainWarp }, neighborDist, seed);
        timing.push({ stage: 'Terrain post-processing', ms: performance.now() - t0 });

        progress(50, 'Deriving plates\u2026');
        t0 = performance.now();
        const { r_plate, plateSeeds, plateIsOcean, plateVec } = deriveSyntheticPlates(mesh, r_elevation);
        timing.push({ stage: 'Synthetic plates', ms: performance.now() - t0 });

        // Classify regions
        const mountain_r = new Set();
        const coastline_r = new Set();
        const ocean_r = new Set();
        for (let r = 0; r < mesh.numRegions; r++) {
            if (r_elevation[r] <= 0) {
                ocean_r.add(r);
            } else if (r_elevation[r] > 0.5) {
                mountain_r.add(r);
            }
            // Coastline: land cell adjacent to ocean
            if (r_elevation[r] > 0) {
                const end = mesh.adjOffset[r + 1];
                for (let ni = mesh.adjOffset[r]; ni < end; ni++) {
                    if (r_elevation[mesh.adjList[ni]] <= 0) {
                        coastline_r.add(r);
                        break;
                    }
                }
            }
        }

        const r_stress = new Float32Array(mesh.numRegions); // no stress for imports
        const debugLayers = { erosionDelta: dl_erosionDelta };
        const nMag = 0;

        let windResult = null, oceanResult = null, precipResult = null, tempResult = null;

        if (!skipClimate) {
            const noise = new SimplexNoise(seed);

            progress(60, 'Simulating wind patterns\u2026');
            t0 = performance.now();
            windResult = computeWind(mesh, r_xyz, r_elevation, plateIsOcean, r_plate, noise);
            timing.push({ stage: 'Wind simulation', ms: performance.now() - t0 });
            debugLayers.pressureSummer = windResult.r_pressure_summer;
            debugLayers.pressureWinter = windResult.r_pressure_winter;
            debugLayers.windSpeedSummer = windResult.r_wind_speed_summer;
            debugLayers.windSpeedWinter = windResult.r_wind_speed_winter;
            debugLayers.continentality = windResult.r_continentality;

            progress(72, 'Computing ocean currents\u2026');
            t0 = performance.now();
            oceanResult = computeOceanCurrents(mesh, r_xyz, r_elevation, windResult);
            timing.push({ stage: 'Ocean currents', ms: performance.now() - t0 });

            progress(80, 'Computing precipitation\u2026');
            t0 = performance.now();
            precipResult = computePrecipitation(mesh, r_xyz, r_elevation, windResult, oceanResult, precipitationOffset, landCoverage);
            timing.push({ stage: 'Precipitation', ms: performance.now() - t0 });
            debugLayers.precipSummer = precipResult.r_precip_summer;
            debugLayers.precipWinter = precipResult.r_precip_winter;
            debugLayers.rainShadowSummer = precipResult.r_rainshadow_summer;
            debugLayers.rainShadowWinter = precipResult.r_rainshadow_winter;

            progress(88, 'Computing temperature\u2026');
            t0 = performance.now();
            tempResult = computeTemperature(mesh, r_xyz, r_elevation, windResult, oceanResult, precipResult, temperatureOffset);
            timing.push({ stage: 'Temperature', ms: performance.now() - t0 });
            debugLayers.tempSummer = tempResult.r_temperature_summer;
            debugLayers.tempWinter = tempResult.r_temperature_winter;
            debugLayers.tempContinentality = tempResult.r_tempContinentality;

            t0 = performance.now();
            debugLayers.koppen = classifyKoppen(mesh, r_elevation, tempResult, precipResult);
            timing.push({ stage: 'Köppen classification', ms: performance.now() - t0 });
        }

        progress(skipClimate ? 75 : 92, 'Computing triangle elevations\u2026');
        t0 = performance.now();
        const t_elevation = computeTriangleElevations(mesh, r_elevation);
        timing.push({ stage: 'Triangle elevations', ms: performance.now() - t0 });

        // Retain state for reapply
        t0 = performance.now();
        W = {
            mesh, r_xyz: new Float32Array(r_xyz), t_xyz: new Float32Array(t_xyz),
            neighborDist,
            r_plate: new Int32Array(r_plate), plateSeeds: new Set(plateSeeds), plateVec,
            plateIsOcean: new Set(plateIsOcean), originalPlateIsOcean: new Set(plateIsOcean),
            plateDensity: {}, plateDensityLand: {}, plateDensityOcean: {},
            prePostElev: new Float32Array(prePostElev),
            r_elevation_final: new Float32Array(r_elevation),
            seed, nMag, noise: new SimplexNoise(seed),
            mountain_r: new Set(mountain_r), coastline_r: new Set(coastline_r), ocean_r: new Set(ocean_r),
            r_stress: new Float32Array(r_stress),
            cachedWind: windResult, cachedOcean: oceanResult
        };
        timing.push({ stage: 'Clone state for retention', ms: performance.now() - t0 });

        const tWorkerTotal = performance.now() - tTotal0;

        // Build result — same shape as handleGenerate's 'done' message
        const result = {
            type: 'done',
            triangles: mesh.triangles,
            halfedges: mesh.halfedges,
            numRegions: mesh.numRegions,
            r_xyz, t_xyz, r_plate,
            plateSeeds: Array.from(plateSeeds),
            plateVec,
            plateIsOcean: Array.from(plateIsOcean),
            originalPlateIsOcean: Array.from(plateIsOcean),
            plateDensity: {}, plateDensityLand: {}, plateDensityOcean: {},
            prePostElev,
            r_elevation, t_elevation,
            mountain_r: Array.from(mountain_r),
            coastline_r: Array.from(coastline_r),
            ocean_r: Array.from(ocean_r),
            r_stress,
            ...buildClimateFields(windResult, oceanResult, precipResult, tempResult),
            skipClimate: !!skipClimate,
            seed, nMag,
            debugLayers,
            _timing: [],
            _pipelineTiming: timing,
            _postTiming: postTiming,
            _workerTotal: tWorkerTotal,
            _params: { N, P: 0, jitter, nMag, numContinents: 0, smoothing, terrainWarp, hydraulicErosion, thermalErosion, ridgeSharpening, glacialErosion, seed }
        };

        const transferList = [
            r_xyz.buffer, t_xyz.buffer, r_plate.buffer,
            prePostElev.buffer, r_elevation.buffer, t_elevation.buffer,
            r_stress.buffer
        ];

        self.postMessage(result, transferList);

    } catch (err) {
        self.postMessage({ type: 'error', message: err.message, stack: err.stack });
    }
}

self.onmessage = (e) => {
    const { cmd } = e.data;
    switch (cmd) {
        case 'generate': handleGenerate(e.data); break;
        case 'reapply': handleReapply(e.data); break;
        case 'editRecompute': handleEditRecompute(e.data); break;
        case 'computeClimate': handleComputeClimate(e.data); break;
        case 'importHeightmap': handleImportHeightmap(e.data); break;
        default: self.postMessage({ type: 'error', message: `Unknown command: ${cmd}` });
    }
};
