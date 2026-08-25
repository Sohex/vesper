#!/usr/bin/env node --test
/**
 * Regression tests for where glacial erosion puts its ice.
 *
 *   node --test tools/test-glacial.mjs
 *
 * Three promises are under test.
 *
 * First, the CONTROL: with no ice mask supplied the placement is bit-identical
 * to the latitude-and-elevation ramp the generator has always used. The
 * reference formula is written out here rather than imported, so the test fails
 * if the ramp is edited and passes only if the two agree number for number.
 *
 * Second, that a supplied mask actually REPLACES the ramp rather than modifying
 * it: an equatorial mask must glaciate the equator, which the ramp never does at
 * any strength, and must be refused outright when its mesh identity disagrees
 * with the mesh it is being indexed into.
 *
 * Third, that the ramp's altitude gate is gravity-INVARIANT on purpose. Both the
 * relief ceiling and a dry-adiabatic freezing height go as 1/g, so a
 * dimensionless gate is the form in which those two cancel. The test pins that,
 * so a later gravity term on the gate has to argue with the algebra in
 * js/glacial-ice.js instead of being added quietly.
 *
 * Needs: npm i delaunator
 */

import test from 'node:test';
import assert from 'node:assert/strict';
import Delaunator from 'delaunator';

import { setDelaunator } from '../js/sphere-mesh.js';
import { buildGlacIdx, validateIceMask } from '../js/glacial-ice.js';
import { runGeneratePipeline } from '../js/pipeline.js';
import { erodeComposite } from '../js/terrain-post.js';
import { computeNeighborDist } from '../js/sphere-mesh.js';
import {
    GLACIAL_LAT_DIVISOR, GLACIAL_ELEV_LOW, GLACIAL_ELEV_HIGH,
    GLACIAL_ELEV_FACTOR_SCALE, GLACIAL_ELEV_FACTOR_LAT_BASE,
    GLACIAL_ELEV_FACTOR_LAT_SCALE,
} from '../js/terrain-config.js';
import { scaledHeightKm } from '../js/color-map.js';
import { hashTypedArray } from '../js/sha256.js';

setDelaunator(Delaunator);

const realLog = console.log;
const quiet = (fn) => { console.log = () => {}; try { return fn(); } finally { console.log = realLog; } };

const PARAMS = {
    N: 8000, P: 24, jitter: 0.75, nMag: 1, numContinents: 5, smoothing: 0,
    hydraulicErosion: 0.5, thermalErosion: 0.3, ridgeSharpening: 0.3,
    glacialErosion: 0.8, terrainWarp: 0, seed: 9091,
    toggledIndices: [], skipClimate: true, lithology: false,
};

let _base = null;
const base = () => (_base ??= quiet(() => runGeneratePipeline(PARAMS)));

/** The placement formula as it stood before glacial-ice.js existed. */
function referenceGlacIdx(mesh, r_xyz, r_elevation, r_isOcean, glacialStrength) {
    const smoothstep = (x, e0, e1) => {
        const t = Math.max(0, Math.min(1, (x - e0) / (e1 - e0)));
        return t * t * (3 - 2 * t);
    };
    const N = mesh.numRegions;
    const idx = new Float32Array(N);
    const thresholdLat = Math.PI / 2 - glacialStrength * Math.PI / GLACIAL_LAT_DIVISOR;
    for (let r = 0; r < N; r++) {
        if (r_isOcean[r]) continue;
        const y = r_xyz[3 * r + 1];
        const polarDist = Math.abs(Math.asin(Math.max(-1, Math.min(1, y))));
        const latFactor = smoothstep(polarDist, thresholdLat, Math.PI / 2);
        const elevFactor = smoothstep(r_elevation[r], GLACIAL_ELEV_LOW, GLACIAL_ELEV_HIGH);
        const latScale = smoothstep(polarDist, Math.PI / 8, Math.PI / 3);
        idx[r] = Math.max(latFactor, elevFactor * GLACIAL_ELEV_FACTOR_SCALE
            * (GLACIAL_ELEV_FACTOR_LAT_BASE + GLACIAL_ELEV_FACTOR_LAT_SCALE * latScale))
            * glacialStrength;
    }
    return idx;
}

function oceanMask(elev) {
    const m = new Uint8Array(elev.length);
    for (let r = 0; r < elev.length; r++) if (elev[r] <= 0) m[r] = 1;
    return m;
}

function latDeg(r_xyz, r) {
    return Math.asin(Math.max(-1, Math.min(1, r_xyz[3 * r + 1]))) * 180 / Math.PI;
}

// ─── Control: the ramp is unchanged ──────────────────────────────────────

test('with no mask the placement is bit-identical to the pre-existing ramp', () => {
    const c = base();
    const isOcean = oceanMask(c.prePostElev);
    for (const strength of [0.3, 0.8, 1.0]) {
        const got = buildGlacIdx(c.mesh, c.r_xyz, c.prePostElev, isOcean, strength).idx;
        const want = referenceGlacIdx(c.mesh, c.r_xyz, c.prePostElev, isOcean, strength);
        assert.equal(hashTypedArray(got), hashTypedArray(want),
            `strength ${strength}: the ramp must not have moved`);
    }
});

test('with no mask the source is reported as the heuristic, not as a measurement', () => {
    const c = base();
    const isOcean = oceanMask(c.prePostElev);
    assert.equal(buildGlacIdx(c.mesh, c.r_xyz, c.prePostElev, isOcean, 0.8).source, 'heuristic');
});

// ─── The gate is dimensionless because the two 1/g cancel ────────────────

test('the ramp is gravity-invariant: reliefScale does not move the altitude gate', () => {
    const c = base();
    const isOcean = oceanMask(c.prePostElev);
    const earth = buildGlacIdx(c.mesh, c.r_xyz, c.prePostElev, isOcean, 0.8,
        { reliefScale: 1 }).idx;
    const heavy = buildGlacIdx(c.mesh, c.r_xyz, c.prePostElev, isOcean, 0.8,
        { reliefScale: 9.80665 / 12.81 }).idx;
    assert.equal(hashTypedArray(earth), hashTypedArray(heavy),
        'a dimensionless gate is where the relief ceiling and the dry-adiabatic '
        + 'freezing height cancel; see the header of js/glacial-ice.js before changing this');
});

test('the gate sits at a physical altitude that scales as 1/g, which is the point', () => {
    // The gate is a fixed model elevation, so the altitude it stands for scales
    // with reliefScale = g_ref/g. A dry-adiabatic freezing height scales the
    // same way. Pinning the first half here means a change that re-anchors the
    // gate to a fixed number of kilometres shows up as a failure.
    const earthKm = scaledHeightKm(GLACIAL_ELEV_LOW, 1, true);
    assert.ok(earthKm > 0);
    for (const g of [4.9, 9.80665, 12.81, 19.6]) {
        const s = 9.80665 / g;
        const km = scaledHeightKm(GLACIAL_ELEV_LOW, s, true);
        assert.ok(Math.abs(km - earthKm * 9.80665 / g) < 1e-9,
            `the gate altitude must go as 1/g; at g=${g} it read ${km}`);
    }
});

// ─── The mask replaces the ramp ──────────────────────────────────────────

test('a supplied mask places the ice, and latitude plays no part', () => {
    const c = base();
    const isOcean = oceanMask(c.prePostElev);
    const N = c.mesh.numRegions;
    // Ice on the equator and nowhere else, which the ramp never produces.
    const values = new Float32Array(N);
    let equatorLand = 0;
    for (let r = 0; r < N; r++) {
        if (isOcean[r]) continue;
        if (Math.abs(latDeg(c.r_xyz, r)) < 10) { values[r] = 1; equatorLand++; }
    }
    assert.ok(equatorLand > 0, 'fixture needs equatorial land');

    const built = buildGlacIdx(c.mesh, c.r_xyz, c.prePostElev, isOcean, 0.8,
        { iceMask: { values, numRegions: N, seed: PARAMS.seed }, seed: PARAMS.seed });
    assert.equal(built.source, 'mask');

    const ramp = referenceGlacIdx(c.mesh, c.r_xyz, c.prePostElev, isOcean, 0.8);
    let maskedEquator = 0, rampEquator = 0, maskedPolar = 0, rampPolar = 0;
    for (let r = 0; r < N; r++) {
        if (isOcean[r]) continue;
        const polar = Math.abs(latDeg(c.r_xyz, r)) > 70;
        if (Math.abs(latDeg(c.r_xyz, r)) < 10) {
            if (built.idx[r] > 0) maskedEquator++;
            if (ramp[r] > 0) rampEquator++;
        }
        if (polar) {
            if (built.idx[r] > 0) maskedPolar++;
            if (ramp[r] > 0) rampPolar++;
        }
    }
    assert.equal(maskedEquator, equatorLand, 'every masked equatorial cell must be glaciated');
    assert.equal(maskedPolar, 0, 'the mask says no polar ice, so there must be none');
    assert.ok(rampPolar > 0, 'the ramp glaciates the poles, which is the difference being tested');
    assert.ok(rampEquator < equatorLand,
        'the ramp does not glaciate the equator, which is why a mask is needed');
});

test('the mask is scaled by the declared strength and never applied to ocean', () => {
    const c = base();
    const isOcean = oceanMask(c.prePostElev);
    const N = c.mesh.numRegions;
    const values = new Float32Array(N).fill(1);
    const idx = buildGlacIdx(c.mesh, c.r_xyz, c.prePostElev, isOcean, 0.4,
        { iceMask: { values, numRegions: N, seed: PARAMS.seed }, seed: PARAMS.seed }).idx;
    for (let r = 0; r < N; r++) {
        assert.equal(idx[r], isOcean[r] ? 0 : Math.fround(0.4),
            `region ${r} must carry the declared strength on land and nothing at sea`);
    }
});

test('a zero mask carves nothing: the terrain matches ice being switched off', () => {
    const c = base();
    const N = c.mesh.numRegions;
    const isOcean = oceanMask(c.prePostElev);
    const neighborDist = computeNeighborDist(c.mesh, c.r_xyz);
    const zeroMask = { values: new Float32Array(N), numRegions: N, seed: PARAMS.seed };

    // Same iteration count on both arms, so the mid-loop priority flood -- which
    // is not glacial and runs for any erosion -- lands identically and the only
    // difference under test is the ice.
    const withZero = Float32Array.from(c.prePostElev);
    const summary = erodeComposite(c.mesh, withZero, c.r_xyz, isOcean,
        0, 0, 0.5, 1.0, 0, 1.0, 0, 8, 0.8, neighborDist, null, null, null,
        { iceMask: zeroMask, seed: PARAMS.seed });
    assert.equal(summary.glacIdxSource, 'mask');
    assert.equal(summary.glaciatedCells, 0);

    const noIce = Float32Array.from(c.prePostElev);
    erodeComposite(c.mesh, noIce, c.r_xyz, isOcean,
        0, 0, 0.5, 1.0, 0, 1.0, 0, 8, 0, neighborDist, null, null, null, null);

    assert.equal(hashTypedArray(withZero), hashTypedArray(noIce),
        'a mask with no ice in it must carve exactly nothing');
});

// ─── The mask is matched by index, so its mesh identity is checked ───────

test('a mask built on another mesh is refused, not silently reindexed', () => {
    const N = 1000;
    const good = { values: new Float32Array(N), numRegions: N, seed: 7 };
    assert.doesNotThrow(() => validateIceMask(good, N, 7));

    assert.throws(() => validateIceMask({ values: new Float32Array(N - 1), seed: 7 }, N, 7),
        /region index/, 'a length mismatch means a different mesh');
    assert.throws(() => validateIceMask({ values: new Float32Array(N), numRegions: N, seed: 8 }, N, 7),
        /seed/, 'a different seed is a different mesh at the same region count');
    assert.throws(() => validateIceMask({ values: new Float32Array(N), numRegions: N + 1, seed: 7 }, N, 7),
        /declares/);
});

test('an out-of-range ice index is refused', () => {
    const N = 8;
    const v = new Float32Array(N);
    v[3] = 1.5;
    assert.throws(() => validateIceMask({ values: v, numRegions: N, seed: 1 }, N, 1), /\[0, 1\]/);
    v[3] = -0.1;
    assert.throws(() => validateIceMask({ values: v, numRegions: N, seed: 1 }, N, 1), /\[0, 1\]/);
    v[3] = NaN;
    assert.throws(() => validateIceMask({ values: v, numRegions: N, seed: 1 }, N, 1), /\[0, 1\]/);
});

// ─── End to end: the mask changes the terrain and says so ────────────────

test('a build carved by a mask is distinguishable from one carved by the ramp', () => {
    const c = base();
    const N = c.mesh.numRegions;
    const isOcean = oceanMask(c.prePostElev);
    const values = new Float32Array(N);
    for (let r = 0; r < N; r++) {
        if (!isOcean[r] && Math.abs(latDeg(c.r_xyz, r)) < 25) values[r] = 1;
    }
    const masked = quiet(() => runGeneratePipeline({
        ...PARAMS,
        iceMask: { values, numRegions: N, seed: PARAMS.seed, provenance: { criterion: 'fixture' } },
    }));

    assert.equal(masked.params.glacialPlacement, 'mask');
    assert.equal(c.params.glacialPlacement, 'heuristic');
    assert.deepEqual(masked.params.iceMask, { criterion: 'fixture' });
    assert.equal(c.params.iceMask, null);
    assert.notEqual(hashTypedArray(masked.r_elevation), hashTypedArray(c.r_elevation),
        'ice placed somewhere else must carve somewhere else');
});
