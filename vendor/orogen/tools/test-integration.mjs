#!/usr/bin/env node --test
/**
 * Regression tests for the pipeline-integration layer: planetary parameters,
 * Gaussian grids, sub-grid orography, hydrology, and planet-code round-tripping.
 *
 *   node --test tools/test-integration.mjs
 *
 * The load-bearing promise here is that **Earth defaults change nothing**. Every
 * one of these features is additive, and a run that names no planet must be
 * bit-identical to one from before planetary parameters existed.
 *
 * Needs: npm i delaunator
 */

import test from 'node:test';
import assert from 'node:assert/strict';
import Delaunator from 'delaunator';

import { setDelaunator } from '../js/sphere-mesh.js';
import { regionLatLon } from '../js/geometry.js';
import {
    gaussianLatitudes, spectralGrid, latitudeEdges, rowForLatitude,
    avgEdgeKm, avgCellAreaKm2, SPECTRAL_GRIDS,
} from '../js/geometry.js';
import { makePlanet, gravityFromMassRadius, EARTH, planetSummary } from '../js/planet-params.js';
import { makeGrid, buildGridLookup, computeSubgridOrography, buildExportBundle } from '../js/data-export.js';
import { computeRivers, computeOceanBasins } from '../js/basins.js';
import { encodePlanetCode, decodePlanetCode } from '../js/planet-code.js';
import { runGeneratePipeline } from '../js/pipeline.js';
import { hashTypedArray } from '../js/sha256.js';
import { collectRegionFields } from '../js/data-export.js';

setDelaunator(Delaunator);

const realLog = console.log;
const quiet = (fn) => { console.log = () => {}; try { return fn(); } finally { console.log = realLog; } };

const PARAMS = {
    N: 20000, P: 24, jitter: 0.75, nMag: 1, numContinents: 5, smoothing: 0,
    hydraulicErosion: 0.5, thermalErosion: 0.3, ridgeSharpening: 0.3,
    glacialErosion: 0.3, terrainWarp: 0, seed: 4242,
    toggledIndices: [], skipClimate: true,
};
let _ctx = null;
const ctx = () => (_ctx ??= quiet(() => runGeneratePipeline(PARAMS)));

// ── Planetary parameters ─────────────────────────────────────────────────

test('THE INVARIANT: Earth defaults change nothing', () => {
    const a = quiet(() => runGeneratePipeline(PARAMS));
    const b = quiet(() => runGeneratePipeline({ ...PARAMS, planet: {} }));
    const c = quiet(() => runGeneratePipeline({ ...PARAMS, planet: { name: 'Earth' } }));
    assert.equal(hashTypedArray(a.r_elevation), hashTypedArray(b.r_elevation));
    assert.equal(hashTypedArray(b.r_elevation), hashTypedArray(c.r_elevation));
    assert.equal(makePlanet({}).reliefScale, 1, 'Earth gravity must give exactly unit relief scale');
});

test('planet validation rejects impossible worlds', () => {
    assert.throws(() => makePlanet({ radiusKm: 0 }), /radius/);
    assert.throws(() => makePlanet({ radiusKm: -1 }), /radius/);
    assert.throws(() => makePlanet({ gravityMS2: 0 }), /gravity/);
    assert.throws(() => makePlanet({ obliquityDeg: 200 }), /obliquity/);
    assert.throws(() => makePlanet({ eccentricity: 1.5 }), /eccentricity/);
    assert.doesNotThrow(() => makePlanet({ radiusKm: 3390, gravityMS2: 3.72 }));   // Mars
});

test('radius rescales every length, gravity rescales relief', () => {
    const p = makePlanet({ radiusKm: 2 * EARTH.radiusKm });
    assert.ok(Math.abs(avgEdgeKm(10000, p.radiusKm) / avgEdgeKm(10000) - 2) < 1e-12);
    assert.ok(Math.abs(avgCellAreaKm2(10000, p.radiusKm) / avgCellAreaKm2(10000) - 4) < 1e-12);

    const heavy = makePlanet({ gravityMS2: 2 * EARTH.gravityMS2 });
    assert.ok(Math.abs(heavy.reliefScale - 0.5) < 1e-12, 'double gravity should halve relief');
});

test('gravity scales PHYSICAL relief as 1/g, leaving the model terrain alone', () => {
    const earth = ctx();
    const heavy = quiet(() => runGeneratePipeline({
        ...PARAMS, planet: { gravityMS2: 2 * EARTH.gravityMS2 },
    }));
    // The model's dimensionless terrain must be untouched: erosion constants are
    // calibrated in model units, and pushing a scaled parameter through the
    // nonlinear height curve damps the effect instead of delivering it.
    assert.equal(hashTypedArray(earth.r_elevation), hashTypedArray(heavy.r_elevation),
        'gravity must not perturb the model terrain');

    const maxKm = (c) => {
        const f = collectRegionFields(c).find(x => x.name === 'elevation_km');
        let m = -Infinity;
        for (const v of f.values) if (v > m) m = v;
        return m;
    };
    const ratio = maxKm(heavy) / maxKm(earth);
    assert.ok(Math.abs(ratio - 0.5) < 1e-6,
        `double gravity gave a height ratio of ${ratio}, expected exactly 0.5`);
});

test('radius rescales exported cell areas to the real sphere', () => {
    const earth = ctx();
    const big = quiet(() => runGeneratePipeline({ ...PARAMS, planet: { radiusKm: 7645.2 } }));
    const totalArea = (c) => {
        const f = collectRegionFields(c).find(x => x.name === 'cell_area');
        let s = 0;
        for (const v of f.values) s += v;
        return s;
    };
    // Within the ~0.1% the centroid-cornered dual cells approximate the sphere.
    for (const [c, r] of [[earth, 6371], [big, 7645.2]]) {
        const expected = 4 * Math.PI * r * r;
        const got = totalArea(c);
        assert.ok(Math.abs(got / expected - 1) < 0.002,
            `exported area ${got} vs sphere ${expected} for R=${r}`);
    }
});

test('elevation and elevation_km are different things, and both are exported', () => {
    const c = ctx();
    const fields = collectRegionFields(c);
    const model = fields.find(f => f.name === 'elevation');
    const km = fields.find(f => f.name === 'elevation_km');
    assert.ok(model && km, 'both elevation fields must be present');
    assert.notEqual(model.units, 'km', 'the model parameter must not claim to be kilometres');
    assert.equal(km.units, 'km');
    // The land/ocean split must agree between them whatever the curve does.
    for (let r = 0; r < c.mesh.numRegions; r++) {
        assert.equal(model.values[r] > 0, km.values[r] > 0, `land/ocean disagree at ${r}`);
    }
    let maxKm = -Infinity;
    for (const v of km.values) if (v > maxKm) maxKm = v;
    assert.ok(maxKm > 3, `max physical height ${maxKm} km looks like the model parameter, not km`);
});

test('gravity from mass and radius matches Earth at Earth values', () => {
    assert.ok(Math.abs(gravityFromMassRadius(1, 1) - EARTH.gravityMS2) < 1e-12);
    assert.ok(gravityFromMassRadius(2, 1) > EARTH.gravityMS2);
    assert.ok(gravityFromMassRadius(1, 2) < EARTH.gravityMS2);
});

test('the summary separates what Orogen uses from what it merely carries', () => {
    const s = planetSummary(makePlanet({ rotationPeriodHours: 12 }));
    assert.deepEqual(s.usedByOrogen, ['radiusKm', 'gravityMS2']);
    assert.ok(s.passedThroughForDownstream.includes('rotationPeriodHours'));
    assert.equal(s.rotationPeriodHours, 12);
});

// ── Erosion scale invariance ─────────────────────────────────────────────

test('THE INVARIANT: erosion strength does not scale with region count', () => {
    // Flow accumulation used to be a cell COUNT and cell spacing a unit-sphere
    // chord, which together made erosive power grow linearly with region count:
    // the same planet at 2.5M regions eroded ~5x harder than at 500k, and basin
    // rims that survived at low resolution were destroyed at high. Both are now
    // normalised to real area and real distance.
    const at = (N) => {
        const c = quiet(() => runGeneratePipeline({ ...PARAMS, N }));
        let land = 0, eroded = 0, maxElev = -Infinity;
        for (let r = 0; r < c.mesh.numRegions; r++) {
            if (c.r_elevation[r] > 0) { land += c.r_elevation[r] * c.cellArea[r]; }
            if (c.r_elevation[r] > maxElev) maxElev = c.r_elevation[r];
            const d = c.debugLayers.erosionDelta[r];
            if (d < 0) eroded += -d * c.cellArea[r];
        }
        return { land, eroded, maxElev };
    };
    const lo = at(20000), hi = at(80000);            // 4x the regions

    const erodedRatio = hi.eroded / lo.eroded;
    assert.ok(erodedRatio > 0.5 && erodedRatio < 2.0,
        `eroded volume changed by ${erodedRatio.toFixed(2)}x for a 4x region count — `
        + 'erosion is scaling with discretisation, not with the terrain');

    const reliefRatio = hi.maxElev / lo.maxElev;
    assert.ok(reliefRatio > 0.8 && reliefRatio < 1.25,
        `peak relief changed by ${reliefRatio.toFixed(2)}x for a 4x region count`);
});

test('the reference configuration keeps its slider calibration', () => {
    // Mean flow seed must be 1 at Earth and the reference region count, or the
    // erosion sliders silently mean something different than they were tuned to.
    const REF_N = 100000, REF_R = 6371;
    const c = quiet(() => runGeneratePipeline({ ...PARAMS, N: REF_N }));
    const refCell = 4 * Math.PI * REF_R * REF_R / REF_N;
    let sum = 0;
    for (let r = 0; r < c.mesh.numRegions; r++) sum += c.cellArea[r] / refCell;
    const mean = sum / c.mesh.numRegions;
    assert.ok(Math.abs(mean - 1) < 0.01, `mean flow seed at the reference is ${mean}, expected 1`);
});

// ── Gaussian grids ───────────────────────────────────────────────────────

test('Gaussian latitudes are symmetric and their weights integrate the sphere', () => {
    for (const n of [8, 32, 64, 96]) {
        const { lat, weights } = gaussianLatitudes(n);
        assert.equal(lat.length, n);
        let sum = 0;
        for (let i = 0; i < n; i++) {
            sum += weights[i];
            assert.ok(Math.abs(lat[i] + lat[n - 1 - i]) < 1e-9, `asymmetric at n=${n}`);
            assert.ok(lat[i] > -90 && lat[i] < 90, `latitude out of range at n=${n}`);
            if (i > 0) assert.ok(lat[i] < lat[i - 1], 'latitudes must descend north to south');
        }
        assert.ok(Math.abs(sum - 2) < 1e-12, `weights sum to ${sum}, expected 2 at n=${n}`);
    }
});

test('Gaussian latitudes differ from equal-angle — which is why this exists', () => {
    const { lat } = gaussianLatitudes(64);
    let maxDiff = 0;
    for (let j = 0; j < 64; j++) {
        maxDiff = Math.max(maxDiff, Math.abs(lat[j] - (90 - (j + 0.5) / 64 * 180)));
    }
    assert.ok(maxDiff > 0.3,
        `max difference from equal-angle is only ${maxDiff}° — the grids would be interchangeable`);
});

test('every Gaussian latitude bins back into its own row', () => {
    for (const n of [32, 64, 128]) {
        const { lat } = gaussianLatitudes(n);
        const edges = latitudeEdges(lat);
        for (let i = 0; i < n; i++) {
            assert.equal(rowForLatitude(Math.sin(lat[i] * Math.PI / 180), edges), i,
                `latitude ${lat[i]} did not bin to row ${i}`);
        }
        assert.equal(edges[0], 1);
        assert.equal(edges[n], -1);
    }
});

test('spectral truncations resolve, and unknown names do not', () => {
    assert.deepEqual(spectralGrid('T42'), { width: 128, height: 64, truncation: 'T42' });
    assert.deepEqual(spectralGrid('t21'), { width: 64, height: 32, truncation: 'T21' });
    assert.equal(spectralGrid('T999'), null);
    assert.equal(spectralGrid('1024x512'), null);
    assert.equal(spectralGrid(null), null);
    for (const [name, g] of Object.entries(SPECTRAL_GRIDS)) {
        assert.equal(g.width, g.height * 2, `${name} must be 2:1`);
    }
});

test('a gaussian grid covers the sphere with no gaps or overlaps', () => {
    const c = ctx();
    const grid = makeGrid(64, 32, 'gaussian');
    const lookup = buildGridLookup(c.mesh, c.r_xyz, grid);
    assert.equal(lookup.length, 64 * 32);
    for (let i = 0; i < lookup.length; i++) {
        assert.ok(lookup[i] >= 0 && lookup[i] < c.mesh.numRegions,
            `cell ${i} resolved to region ${lookup[i]}`);
    }
});

// ── Sub-grid orography ───────────────────────────────────────────────────

test('sub-grid statistics bracket the mean and are non-negative where defined', () => {
    const c = ctx();
    const grid = makeGrid(64, 32, 'gaussian');
    const lookup = buildGridLookup(c.mesh, c.r_xyz, grid);
    const { lat, lon } = regionLatLon(c.r_xyz);
    const st = computeSubgridOrography(c.mesh, c.r_xyz, c.r_elevation, grid, lookup,
        lat, lon, c.cellArea, c.neighborDist, c.planet.radiusKm);
    const { orog_mean, orog_std, orog_min, orog_max, orog_count, orog_anisotropy } = st.fields;
    for (let i = 0; i < orog_mean.length; i++) {
        if (orog_count[i] === 0) continue;
        assert.ok(orog_std[i] >= 0, `negative std at ${i}`);
        assert.ok(orog_min[i] <= orog_mean[i] + 1e-5, `min above mean at ${i}`);
        assert.ok(orog_max[i] >= orog_mean[i] - 1e-5, `max below mean at ${i}`);
        assert.ok(orog_anisotropy[i] >= 0 && orog_anisotropy[i] <= 1,
            `anisotropy ${orog_anisotropy[i]} out of range at ${i}`);
    }
});

test('orog_* is physical km, not the model shaping parameter', () => {
    // The regression: computeSubgridOrography was wired to data.r_elevation, the
    // dimensionless model field, while declaring units 'km'. orog_mean came out
    // bit-identical to `elevation` and read a 4.56 km peak as 1.14 "km". It is a
    // wiring bug, so it has to be caught through buildExportBundle rather than
    // by calling computeSubgridOrography directly with a correct array.
    const gridOf = (bundle, name) => {
        const f = bundle.files.find(x => x.name === `grid/${name}.bin`);
        assert.ok(f, `${name} missing from the bundle`);
        return new Float32Array(f.bytes.buffer, f.bytes.byteOffset, f.bytes.byteLength / 4);
    };
    const bundle = buildExportBundle(ctx(), {
        raw: false, gridWidth: 32, gridHeight: 16, subgridOrography: true,
    });
    const orogMean = gridOf(bundle, 'orog_mean');
    const elevKm = gridOf(bundle, 'elevation_km');
    const elevModel = gridOf(bundle, 'elevation');
    const count = gridOf(bundle, 'orog_count');

    let checked = 0, matchesModel = 0;
    for (let i = 0; i < orogMean.length; i++) {
        if (count[i] === 0) continue;
        checked++;
        assert.ok(Math.abs(orogMean[i] - elevKm[i]) < 1e-4,
            `orog_mean ${orogMean[i]} is not the km elevation ${elevKm[i]} at cell ${i}`);
        if (Math.abs(orogMean[i] - elevModel[i]) < 1e-6) matchesModel++;
    }
    assert.ok(checked > 50, 'too few populated cells to be a meaningful check');
    // Sea-level cells are 0 in both conventions, so demand only that the two
    // fields are not the SAME field — under the bug every cell matched.
    assert.ok(matchesModel < checked,
        'orog_mean is bit-identical to the model elevation field: the km conversion was skipped');
});

test('orog_* carries the 1/g relief scaling, like elevation_km', () => {
    const bundleAt = (gravityMS2) => buildExportBundle(
        quiet(() => runGeneratePipeline({ ...PARAMS, planet: { gravityMS2 } })),
        { raw: false, gridWidth: 32, gridHeight: 16, subgridOrography: true });
    const maxOf = (bundle) => {
        const f = bundle.files.find(x => x.name === 'grid/orog_max.bin');
        const a = new Float32Array(f.bytes.buffer, f.bytes.byteOffset, f.bytes.byteLength / 4);
        let m = -Infinity;
        for (let i = 0; i < a.length; i++) if (a[i] > m) m = a[i];
        return m;
    };
    const earth = maxOf(bundleAt(EARTH.gravityMS2));
    const heavy = maxOf(bundleAt(2 * EARTH.gravityMS2));
    assert.ok(earth > 0.5, `expected real relief on the Earth-gravity run, got ${earth}`);
    // Land height scales as 1/g, so the highest point halves at double gravity.
    assert.ok(Math.abs(heavy / earth - 0.5) < 1e-5,
        `peak orog_max should halve at 2g: ${earth} -> ${heavy}`);
});

test('flat terrain has zero sub-grid roughness', () => {
    const c = ctx();
    const grid = makeGrid(32, 16, 'uniform');
    const lookup = buildGridLookup(c.mesh, c.r_xyz, grid);
    const { lat, lon } = regionLatLon(c.r_xyz);
    const flat = new Float32Array(c.mesh.numRegions).fill(0.5);
    const st = computeSubgridOrography(c.mesh, c.r_xyz, flat, grid, lookup,
        lat, lon, c.cellArea, c.neighborDist, c.planet.radiusKm);
    for (let i = 0; i < st.fields.orog_std.length; i++) {
        if (st.fields.orog_count[i] === 0) continue;
        assert.ok(st.fields.orog_std[i] < 1e-5, `std ${st.fields.orog_std[i]} on flat ground`);
    }
});

// ── Hydrology ────────────────────────────────────────────────────────────

test('flow accumulation is monotonic downstream and starts at cell area', () => {
    const c = ctx();
    const { flowAccumKm2, drainTo } = c.hydrology
        ? { flowAccumKm2: c.hydrology.r_flowAccumKm2, drainTo: c.hydrology.r_drainTo }
        : computeRivers(c.mesh, c.r_elevation, new Uint8Array(c.mesh.numRegions), c.cellArea);
    for (let r = 0; r < c.mesh.numRegions; r++) {
        assert.ok(flowAccumKm2[r] >= c.cellArea[r] - 1e-3,
            `accumulation ${flowAccumKm2[r]} below own area ${c.cellArea[r]} at ${r}`);
        const t = drainTo[r];
        if (t >= 0) {
            assert.ok(flowAccumKm2[t] >= flowAccumKm2[r] - 1e-3,
                `accumulation decreases downstream at ${r} -> ${t}`);
        }
    }
});

test('river mouths discharge into the open ocean, largest first', () => {
    const c = ctx();
    const H = c.hydrology;
    assert.ok(H.mouthCount > 0, 'no river mouths found');
    for (let i = 1; i < H.mouths.length; i++) {
        assert.ok(H.mouths[i - 1].dischargeAreaKm2 >= H.mouths[i].dischargeAreaKm2,
            'mouths must be sorted by discharge area');
    }
    for (const m of H.mouths) {
        assert.ok(c.r_elevation[m.region] > 0, `mouth ${m.region} is not a land cell`);
        assert.ok(m.dischargeAreaKm2 > 0);
    }
});

test('marginal seas are disconnected from the world ocean and report a sill', () => {
    const c = ctx();
    const isOcean = new Uint8Array(c.mesh.numRegions);
    for (let r = 0; r < c.mesh.numRegions; r++) if (c.r_elevation[r] <= 0) isOcean[r] = 1;
    const { basins, isOpenOcean } = computeOceanBasins(c.mesh, c.r_elevation, isOcean, c.cellArea);
    for (const b of basins) {
        assert.ok(b.areaKm2 > 0);
        assert.ok(b.maxDepthKm >= 0, `negative max depth on ${b.id}`);
        for (const r of b.members) {
            assert.equal(isOpenOcean[r], 0, 'a marginal sea cell must not be open ocean');
            assert.ok(c.r_elevation[r] <= 0, 'a marginal sea cell must be submerged');
        }
    }
});

// ── Planet code ──────────────────────────────────────────────────────────

const CODE_ARGS = [16236323, 510000, 0.75, 100, 10, 0.4, 0.75, 0.1, 0.8, 0.6, 0.3, 0.5, 0.75, 0.85, 0, 0, 0.45];

test('THE REGRESSION: planet codes round-trip the basin and lithology sliders', () => {
    for (const [bas, lit] of [[0, 0], [3, 1.0], [8, 0.5], [5, 0.25], [1, 0.05]]) {
        const code = encodePlanetCode(...CODE_ARGS, bas, lit);
        const d = decodePlanetCode(code);
        assert.ok(d, `failed to decode ${code}`);
        assert.equal(d.basinSlider, bas, `basin slider lost in ${code}`);
        assert.ok(Math.abs(d.lithologyStrength - lit) < 1e-9, `rock contrast lost in ${code}`);
    }
});

test('planet codes round-trip every other slider unchanged', () => {
    const code = encodePlanetCode(...CODE_ARGS, 3, 1.0);
    const d = decodePlanetCode(code);
    assert.equal(d.seed, CODE_ARGS[0]);
    assert.equal(d.N, CODE_ARGS[1]);
    assert.equal(d.P, CODE_ARGS[3]);
    assert.equal(d.numContinents, CODE_ARGS[4]);
    assert.ok(Math.abs(d.roughness - 0.4) < 1e-9);
    assert.ok(Math.abs(d.terrainWarp - 0.75) < 1e-9);
    assert.ok(Math.abs(d.glacialErosion - 0.8) < 1e-9);
    assert.ok(Math.abs(d.soilCreep - 0.75) < 1e-9);
    assert.ok(Math.abs(d.continentSizeVariety - 0.85) < 1e-9);
    assert.ok(Math.abs(d.landCoverage - 0.45) < 1e-9);
});

test('plate toggles survive alongside the new sliders', () => {
    const code = encodePlanetCode(...CODE_ARGS, 3, 1.0, [1, 5, 9]);
    const d = decodePlanetCode(code);
    assert.deepEqual(d.toggledIndices, [1, 5, 9]);
    assert.equal(d.basinSlider, 3);
});

test('legacy codes still decode, to the behaviour those planets actually had', () => {
    const d = decodePlanetCode('09oa7lj8kek17v5639gvhe');   // a real 22-char code
    assert.ok(d, 'legacy code failed to decode');
    assert.equal(d.seed, 16236323);
    assert.equal(d.N, 510000);
    assert.equal(d.P, 100);
    assert.ok(Math.abs(d.landCoverage - 0.45) < 1e-9);
    // Predates both sliders, so it must decode to upstream behaviour rather
    // than to today's defaults — otherwise the planet silently changes.
    assert.equal(d.basinSlider, 0, 'legacy code must not gain basin preservation');
    assert.equal(d.lithologyStrength, 0, 'legacy code must not gain rock contrast');
});

test('malformed codes are rejected rather than mis-decoded', () => {
    for (const bad of ['', 'zzzz', 'not-a-code', '!!!', '0'.repeat(80)]) {
        assert.equal(decodePlanetCode(bad), null, `accepted "${bad}"`);
    }
    assert.equal(decodePlanetCode(null), null);
    assert.equal(decodePlanetCode(42), null);
});

test('encoding an out-of-range value throws instead of silently wrapping', () => {
    assert.throws(() => encodePlanetCode(...CODE_ARGS, 99, 0), /outside its slider range/);
});
