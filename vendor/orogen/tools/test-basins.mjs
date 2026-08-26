#!/usr/bin/env node --test
/**
 * Regression tests for endorheic basin preservation.
 *
 *   node --test tools/test-basins.mjs
 *
 * The headline test is "preserved basins keep their rims": the whole feature is
 * a promise that the drainage conditioning will not carve an outlet through a
 * closed basin, and the previous attempt at this failed silently because
 * nothing checked it. Most of the rest exist to make sure that promise cannot
 * quietly regress.
 *
 * Needs: npm i delaunator
 */

import test from 'node:test';
import assert from 'node:assert/strict';
import Delaunator from 'delaunator';

import { setDelaunator, buildSphere, generateTriangleCenters, computeNeighborDist } from '../js/sphere-mesh.js';
import { makeRng } from '../js/rng.js';
import { EARTH } from '../js/planet-params.js';
import { BASIN_MIN_DEPTH_KM } from '../js/terrain-config.js';
import {
    detectBasins, selectBasins, buildBasinProtection, basinId, attachHypsometry,
    parseBasinList, inciseOutlets,
    computeFillSurface, findOpenOcean, remeasureBasins,
    computeFinalTerrainDrainage, computeFinalPreservedCatchments, compareCatchments,
    applyInlandWaterLevels, surfaceClasses, basinHypsometry,
    BASIN_SURFACE_CLASSES, usesLandHeightBranch,
    auditCatchmentConsistency, computeRivers,
} from '../js/basins.js';
import { erodeComposite, priorityFloodCarve, assertDividesNotLowered } from '../js/terrain-post.js';
import { regionCellArea, uniformCellArea } from '../js/geometry.js';
import { runGeneratePipeline } from '../js/pipeline.js';
import { sha256, hashTypedArray, hashJson } from '../js/sha256.js';

setDelaunator(Delaunator);

// Silence the simulation's own diagnostics.
const realLog = console.log;
const quiet = (fn) => { console.log = () => {}; try { return fn(); } finally { console.log = realLog; } };

// ── Shared fixtures ──────────────────────────────────────────────────────

const PARAMS = {
    N: 20000, P: 24, jitter: 0.75, nMag: 1, numContinents: 5, smoothing: 0,
    hydraulicErosion: 0.4, thermalErosion: 0.3, ridgeSharpening: 0.3,
    glacialErosion: 0.3, terrainWarp: 0, seed: 4242,
    toggledIndices: [], skipClimate: true,
};

let _ctx = null;
const ctx = () => (_ctx ??= quiet(() => runGeneratePipeline({ ...PARAMS, preserveBasins: true })));

/** A small synthetic world: one crater basin on a domed island. */
function craterWorld() {
    const { mesh, r_xyz } = buildSphere(4000, 0.75, makeRng(7));
    const n = mesh.numRegions;
    const elev = new Float32Array(n);
    // Island centred on +y, with a circular moat (the basin) inside a rim.
    for (let r = 0; r < n; r++) {
        const y = r_xyz[3 * r + 1];
        const d = Math.acos(Math.max(-1, Math.min(1, y)));   // angular distance from the pole
        if (d > 0.9) { elev[r] = -1; continue; }             // ocean
        if (d > 0.55) { elev[r] = 0.4; continue; }           // outer flank
        if (d > 0.35) { elev[r] = 1.2; continue; }           // rim
        elev[r] = 0.3;                                       // basin floor, 900 m below the rim
    }
    const isOcean = new Uint8Array(n);
    for (let r = 0; r < n; r++) if (elev[r] <= 0) isOcean[r] = 1;
    return { mesh, r_xyz, elev, isOcean, t_xyz: generateTriangleCenters(mesh, r_xyz) };
}

// ── Non-mutating diagnostics ─────────────────────────────────────────────

test('detectBasins does not mutate the terrain it measures', () => {
    const w = craterWorld();
    const before = Float32Array.from(w.elev);
    detectBasins(w.mesh, w.elev, w.isOcean, {});
    assert.deepEqual(Array.from(w.elev), Array.from(before));
});

test('computeFillSurface never lowers a cell below its own elevation', () => {
    const w = craterWorld();
    const { isOpenOcean } = findOpenOcean(w.mesh, w.isOcean);
    const fill = computeFillSurface(w.mesh, w.elev, isOpenOcean);
    for (let r = 0; r < w.mesh.numRegions; r++) {
        assert.ok(fill[r] >= w.elev[r] - 1e-6, `fill[${r}] below elevation`);
    }
});

// ── Detection ────────────────────────────────────────────────────────────

test('detects the synthetic crater with the right depth and spill level', async () => {
    const w = craterWorld();
    const basins = detectBasins(w.mesh, w.elev, w.isOcean, { cellArea: uniformCellArea(w.mesh.numRegions) });
    assert.ok(basins.length > 0, 'no depression found in a world with an obvious crater');
    const crater = basins.reduce((a, b) => (b.volumeKm3 > a.volumeKm3 ? b : a));
    // Unsuffixed fields are the MODEL parameter the terrain is built in.
    // Elevations round-trip through Float32Array, so compare with a tolerance.
    assert.ok(Math.abs(crater.sinkElevation - 0.3) < 1e-6, `sink ${crater.sinkElevation} != 0.3`);
    assert.ok(Math.abs(crater.spillElevation - 1.2) < 1e-6, `spill ${crater.spillElevation} != 1.2`);
    assert.ok(Math.abs(crater.depth - 0.9) < 1e-6, `depth ${crater.depth} != 0.9`);
    // ...and the Km-suffixed ones are physical kilometres, via the LAND branch
    // of the height curve. A field named Km must be km.
    const { elevToHeightKm } = await import('../js/color-map.js');
    assert.ok(Math.abs(crater.sinkElevationKm - elevToHeightKm(0.3, true)) < 1e-5);
    assert.ok(Math.abs(crater.spillElevationKm - elevToHeightKm(1.2, true)) < 1e-5);
    assert.ok(crater.depthKm > crater.depth,
        'the height curve is convex here, so physical depth exceeds the model difference');
    assert.equal(crater.spillsInto, 'ocean');
});

test('an inland sea below sea level is a basin, not part of the ocean', () => {
    const w = craterWorld();
    // Drop the crater floor below sea level; it is still enclosed by the rim.
    for (let r = 0; r < w.mesh.numRegions; r++) if (Math.abs(w.elev[r] - 0.3) < 1e-6) w.elev[r] = -0.2;
    const isOcean = new Uint8Array(w.mesh.numRegions);
    for (let r = 0; r < w.mesh.numRegions; r++) if (w.elev[r] <= 0) isOcean[r] = 1;
    const { isOpenOcean } = findOpenOcean(w.mesh, isOcean);
    const floorCells = [];
    for (let r = 0; r < w.mesh.numRegions; r++) if (Math.abs(w.elev[r] + 0.2) < 1e-6) floorCells.push(r);
    assert.ok(floorCells.length > 0);
    for (const r of floorCells) {
        assert.equal(isOpenOcean[r], 0, 'enclosed below-sea-level floor must not count as open ocean');
    }
});

test('basin ids are deterministic and depend on both sink and membership', () => {
    assert.equal(basinId(5, [1, 2, 3]), basinId(5, [1, 2, 3]));
    assert.notEqual(basinId(5, [1, 2, 3]), basinId(6, [1, 2, 3]));
    assert.notEqual(basinId(5, [1, 2, 3]), basinId(5, [1, 2, 4]));
    assert.match(basinId(5, [1, 2, 3]), /^orogen-basin-r[0-9a-z]+-[0-9a-z]+$/);
});

// ── Selection ────────────────────────────────────────────────────────────

test('selection floors are physical, and explicit ids bypass them', () => {
    const basins = [
        { id: 'a', index: 0, parentIndex: -1, nestDepth: 0, depth: 1.0, depthKm: 1.0, areaKm2: 50000, cellCount: 500, volumeKm3: 100, sink: 1 },
        { id: 'b', index: 1, parentIndex: -1, nestDepth: 0, depth: 0.001, depthKm: 0.001, areaKm2: 50000, cellCount: 500, volumeKm3: 50, sink: 2 },
        { id: 'c', index: 2, parentIndex: -1, nestDepth: 0, depth: 1.0, depthKm: 1.0, areaKm2: 10, cellCount: 500, volumeKm3: 10, sink: 3 },
        { id: 'd', index: 3, parentIndex: -1, nestDepth: 0, depth: 1.0, depthKm: 1.0, areaKm2: 50000, cellCount: 2, volumeKm3: 5, sink: 4 },
    ];
    const { selected } = selectBasins(basins, {});
    assert.deepEqual(selected.map(b => b.id), ['a'], 'only the basin clearing every floor');

    const forced = selectBasins(basins, { explicitIds: ['c', 'd'] });
    assert.deepEqual(forced.selected.map(b => b.id).sort(), ['a', 'c', 'd']);
    assert.equal(forced.selected.find(b => b.id === 'c').selectedBy, 'explicit');
});

test('the depth floor is compared in kilometres, not in model units', () => {
    // The two currencies disagree wherever the height curve is not the
    // identity, and BASIN_MIN_DEPTH_KM is declared, named and published as a
    // physical depth. A fixture whose `depth` and `depthKm` agree cannot tell
    // which one selectBasins read, so this one makes them disagree in both
    // directions and each half can fail on its own.
    const basins = [
        // Deep in model units, far too shallow to be a landform. On high ground
        // the quartic branch of the curve turns a small model-unit drop into a
        // large physical one; this is the opposite case, and the saturation at
        // model elevation 1 makes the extreme of it, where a depression's
        // physical depth -- and so its lake capacity -- is exactly zero.
        { id: 'model-deep', index: 0, parentIndex: -1, nestDepth: 0, depth: 1.0, depthKm: 0.0, areaKm2: 50000, cellCount: 500, volumeKm3: 100, sink: 1 },
        // Shallow in model units, a real 400 m depression.
        { id: 'km-deep', index: 1, parentIndex: -1, nestDepth: 0, depth: 0.001, depthKm: 0.4, areaKm2: 50000, cellCount: 500, volumeKm3: 50, sink: 2 },
    ];
    assert.deepEqual(selectBasins(basins, {}).selected.map(b => b.id), ['km-deep'],
        'the physical depth decides, so a zero-capacity depression is out and a 400 m one is in');
});

test('the depth floor is compared at REFERENCE gravity, not at the planet\'s', () => {
    // `depthKm` carries the 1/g relief scaling and `selectionDepthKm` does not.
    // The floor is compared against the second, because selectBasins runs during
    // generation and what it returns is notched and protected into r_elevation:
    // a selection that moved with gravity would move the model terrain with it.
    // The two fields disagree in both directions here -- a heavy planet, where
    // the scaled depth is the smaller, and a light one, where it is the larger --
    // so each half fails on its own if the comparison reads `depthKm`.
    const basins = [
        { id: 'heavy-planet', index: 0, parentIndex: -1, nestDepth: 0, depth: 1.0, depthKm: 0.01, selectionDepthKm: 0.4, areaKm2: 50000, cellCount: 500, volumeKm3: 100, sink: 1 },
        { id: 'light-planet', index: 1, parentIndex: -1, nestDepth: 0, depth: 1.0, depthKm: 0.4, selectionDepthKm: 0.01, areaKm2: 50000, cellCount: 500, volumeKm3: 50, sink: 2 },
    ];
    assert.deepEqual(selectBasins(basins, {}).selected.map(b => b.id), ['heavy-planet'],
        'the reference-gravity depth decides, whichever way the scaling runs');
});

test('nested basins are not auto-selected, and an unknown id warns', () => {
    const basins = [
        { id: 'outer', index: 0, parentIndex: -1, nestDepth: 0, depth: 1, depthKm: 1, areaKm2: 50000, cellCount: 500, volumeKm3: 100, sink: 1 },
        { id: 'inner', index: 1, parentIndex: 0, nestDepth: 1, depth: 1, depthKm: 1, areaKm2: 50000, cellCount: 500, volumeKm3: 50, sink: 2 },
    ];
    assert.deepEqual(selectBasins(basins, {}).selected.map(b => b.id), ['outer']);
    assert.equal(selectBasins(basins, { includeNested: true }).selected.length, 2);

    const nestedWarn = selectBasins(basins, { explicitIds: ['inner'] });
    assert.equal(nestedWarn.selected.length, 2);
    assert.ok(nestedWarn.warnings.some(w => /nested inside/.test(w)));

    const unknown = selectBasins(basins, { explicitIds: ['nope'] });
    assert.ok(unknown.warnings.some(w => /No basin with id "nope"/.test(w)));
});

test('selection is disabled wholesale by enabled:false but explicit ids still win', () => {
    const basins = [{ id: 'a', index: 0, parentIndex: -1, nestDepth: 0, depth: 1, depthKm: 1, areaKm2: 50000, cellCount: 500, volumeKm3: 1, sink: 1 }];
    assert.equal(selectBasins(basins, { enabled: false }).selected.length, 0);
    assert.equal(selectBasins(basins, { enabled: false, explicitIds: ['a'] }).selected.length, 1);
});

// ── Protection ───────────────────────────────────────────────────────────

test('protection marks sinks terminal and builds a divide ring around members', () => {
    const w = craterWorld();
    const basins = detectBasins(w.mesh, w.elev, w.isOcean, {});
    const crater = basins.reduce((a, b) => (b.volumeKm3 > a.volumeKm3 ? b : a));
    const p = buildBasinProtection(w.mesh, [crater], { ringCells: 2 });

    assert.equal(p.terminalRoot[crater.sink], 1);
    assert.equal(p.noLower[crater.sink], 0, 'the sink itself must stay erodible');
    assert.ok(p.noLower.reduce((a, b) => a + b, 0) > 0, 'divide ring is empty');
    for (const r of crater.members) assert.equal(p.basinIndex[r], 0);
});

test('THE INVARIANT: the guard throws when a divide is lowered, and not otherwise', () => {
    const before = new Float32Array([1.0, 2.0, 3.0]);
    const same = new Float32Array([1.0, 2.0, 3.0]);
    const raised = new Float32Array([1.5, 2.0, 3.0]);
    const lowered = new Float32Array([1.0, 1.5, 3.0]);

    assert.doesNotThrow(() => assertDividesNotLowered([0, 1, 2], before, same));
    assert.doesNotThrow(() => assertDividesNotLowered([0, 1, 2], before, raised),
        'raising a divide is allowed — Pass 3 only ever raises');
    assert.throws(() => assertDividesNotLowered([0, 1, 2], before, lowered),
        /lowered protected divide cell 1/);
    // Below the epsilon is not a breach.
    assert.doesNotThrow(() => assertDividesNotLowered([0], before, new Float32Array([1.0 - 1e-9, 2, 3])));
});

test('priorityFloodCarve does not lower any protected divide, even at max carve', () => {
    const w = craterWorld();
    const basins = detectBasins(w.mesh, w.elev, w.isOcean, {});
    const crater = basins.reduce((a, b) => (b.volumeKm3 > a.volumeKm3 ? b : a));
    const p = buildBasinProtection(w.mesh, [crater]);

    // Deepen the basin so the carve has a large deficit it would love to
    // redistribute straight through the rim.
    const elev = Float32Array.from(w.elev);
    for (const r of crater.members) elev[r] -= 0.5;
    const divides = [];
    for (let r = 0; r < w.mesh.numRegions; r++) if (p.noLower[r]) divides.push(r);
    const before = Float32Array.from(divides.map(r => elev[r]));
    assert.ok(divides.length > 0, 'no divide cells to test');

    priorityFloodCarve(w.mesh, elev, w.isOcean, 0.85, p);

    for (let i = 0; i < divides.length; i++) {
        assert.ok(elev[divides[i]] >= before[i] - 1e-7,
            `divide cell ${divides[i]} fell from ${before[i]} to ${elev[divides[i]]}`);
    }
});

test('protecting every cell forces filling instead of carving', () => {
    const w = craterWorld();
    const basins = detectBasins(w.mesh, w.elev, w.isOcean, {});
    const crater = basins.reduce((a, b) => (b.volumeKm3 > a.volumeKm3 ? b : a));
    const p = buildBasinProtection(w.mesh, [crater]);
    const all = { terminalRoot: p.terminalRoot, noLower: new Uint8Array(w.mesh.numRegions).fill(1) };
    const elev = Float32Array.from(w.elev);
    for (const r of crater.members) elev[r] -= 0.5;
    const before = Float32Array.from(elev);

    assert.doesNotThrow(() => priorityFloodCarve(w.mesh, elev, w.isOcean, 0.85, all));
    for (let r = 0; r < w.mesh.numRegions; r++) {
        if (w.isOcean[r]) continue;
        assert.ok(elev[r] >= before[r] - 1e-7, `cell ${r} was lowered despite blanket protection`);
    }
});

test('a protected divide survives a full erosion run', () => {
    const w = craterWorld();
    const basins = detectBasins(w.mesh, w.elev, w.isOcean, {});
    const crater = basins.reduce((a, b) => (b.volumeKm3 > a.volumeKm3 ? b : a));
    const p = buildBasinProtection(w.mesh, [crater]);
    const neighborDist = computeNeighborDist(w.mesh, w.r_xyz);

    const withProt = Float32Array.from(w.elev);
    const without = Float32Array.from(w.elev);
    quiet(() => {
        erodeComposite(w.mesh, withProt, w.r_xyz, w.isOcean, 8, 0.00024, 0.5, 1.0, 3, 1.08, 0.045, 3, 0.3, neighborDist, p);
        erodeComposite(w.mesh, without,  w.r_xyz, w.isOcean, 8, 0.00024, 0.5, 1.0, 3, 1.08, 0.045, 3, 0.3, neighborDist, null);
    });

    const cellArea = regionCellArea(w.mesh, w.r_xyz, w.t_xyz);
    const oceanOf = (e) => { const m = new Uint8Array(e.length); for (let r = 0; r < e.length; r++) if (e[r] <= 0) m[r] = 1; return m; };
    const kept = remeasureBasins(w.mesh, withProt, oceanOf(withProt), [crater], cellArea)[0];
    const lost = remeasureBasins(w.mesh, without,  oceanOf(without),  [crater], cellArea)[0];

    assert.ok(kept.retainedFraction > 0.8,
        `protected basin kept only ${(100 * kept.retainedFraction).toFixed(0)}% of its depth`);
    assert.ok(kept.retainedFraction > lost.retainedFraction * 2,
        `protection made little difference: ${kept.retainedFraction} vs ${lost.retainedFraction}`);
});

// ── Hypsometry ───────────────────────────────────────────────────────────

test('hypsometry is monotonic in level, area and volume', () => {
    const w = craterWorld();
    const cellArea = uniformCellArea(w.mesh.numRegions);
    const basins = detectBasins(w.mesh, w.elev, w.isOcean, { cellArea });
    const crater = basins.reduce((a, b) => (b.volumeKm3 > a.volumeKm3 ? b : a));
    // Hypsometry is attached after selection, not during detection — computing
    // it for every depression on a large planet would be almost entirely waste.
    assert.equal(crater.hypsometry, undefined, 'detection must not compute hypsometry');
    attachHypsometry([crater], w.elev, cellArea);
    const h = crater.hypsometry;
    assert.ok(h.length >= 2);
    for (let i = 1; i < h.length; i++) {
        assert.ok(h[i].levelKm >= h[i - 1].levelKm, 'level must ascend');
        assert.ok(h[i].areaKm2 >= h[i - 1].areaKm2 - 1e-6, 'flooded area must not shrink as level rises');
        assert.ok(h[i].volumeKm3 >= h[i - 1].volumeKm3 - 1e-6, 'volume must not shrink as level rises');
    }
    assert.ok(Math.abs(h[0].depthKm) < 1e-9, 'curve starts at the sink');
    assert.ok(Math.abs(h[h.length - 1].levelKm - crater.spillElevationKm) < 1e-6, 'curve ends at the spill point');
});

test('basinHypsometry handles a degenerate flat basin without dividing by zero', () => {
    const elev = new Float32Array([0.5, 0.5, 0.5]);
    const h = basinHypsometry(elev, [0, 1, 2], new Float32Array([1, 1, 1]), 0.5, 0.5, 4);
    assert.equal(h.length, 4);
    for (const p of h) assert.ok(Number.isFinite(p.areaKm2) && Number.isFinite(p.volumeKm3));
});

// ── Water levels ─────────────────────────────────────────────────────────

test('water levels are validated, clamped at the spill point, and never inferred', () => {
    const w = craterWorld();
    const basins = detectBasins(w.mesh, w.elev, w.isOcean, {});
    const crater = basins.reduce((a, b) => (b.volumeKm3 > a.volumeKm3 ? b : a));

    // No level supplied ⇒ dry basin. Orogen must not invent one.
    const dry = applyInlandWaterLevels(w.mesh, w.elev, [crater], {});
    assert.equal(dry.applied[0].levelKm, null);
    assert.equal(dry.isInlandWater.reduce((a, b) => a + b, 0), 0);

    const partial = applyInlandWaterLevels(w.mesh, w.elev, [crater], { [crater.id]: 0.6 });
    assert.equal(partial.applied[0].levelKm, 0.6);
    assert.ok(partial.isInlandWater.reduce((a, b) => a + b, 0) > 0);

    const over = applyInlandWaterLevels(w.mesh, w.elev, [crater], { [crater.id]: 99 });
    assert.equal(over.applied[0].levelKm, crater.spillElevationKm, 'must clamp to the spill point');
    assert.ok(over.warnings.some(x => /above its spill point/.test(x)));

    const under = applyInlandWaterLevels(w.mesh, w.elev, [crater], { [crater.id]: -5 });
    assert.equal(under.applied[0].levelKm, null);
    assert.ok(under.warnings.some(x => /below its floor/.test(x)));

    assert.ok(applyInlandWaterLevels(w.mesh, w.elev, [crater], { nope: 1 })
        .warnings.some(x => /not a preserved basin/.test(x)));
});

test('the land height branch is an explicit set, not the complement of ocean', () => {
    // A set defined as "not ocean" is correct until the vocabulary grows, then
    // the new member joins it silently. Downstream lost a rock class exactly
    // that way. Assert the membership positively, per class.
    assert.equal(usesLandHeightBranch(0), false, 'open ocean is seabed');
    assert.equal(usesLandHeightBranch(1), true, 'subaerial land');
    assert.equal(usesLandHeightBranch(2), true, 'an inland lake floor is continental crust');
    // Any class this test does not know about must be an explicit decision, not
    // an inheritance. If BASIN_SURFACE_CLASSES grows, this fails until someone
    // decides which branch the new class takes.
    assert.deepEqual(Object.keys(BASIN_SURFACE_CLASSES).map(Number).sort(), [0, 1, 2],
        'a surface class was added — decide its height branch, do not let it default');
});

test('the relief curve is strictly increasing over the whole range the model produces', async () => {
    // The shape function t^4(5-4t) is a Hermite interpolant on [0, 1] and turns
    // over above it, so the curve used to clamp at 1 and publish a quarter of a
    // per cent of land at one identical height. A closed basin with both its
    // sink and its spill in that band then had a depth of exactly zero.
    // The model's elevation parameter is not bounded by 1 and reaches past 1.4
    // on the registered builds, so the curve is asserted well beyond that.
    const { elevToHeightKm, RELIEF_KM_PER_UNIT } = await import('../js/color-map.js');

    let prev = -Infinity;
    for (let t = 0.001; t <= 2.0; t += 0.001) {
        const h = elevToHeightKm(t, true);
        assert.ok(h > prev, `height must increase at elev ${t.toFixed(3)}: ${h} <= ${prev}`);
        prev = h;
    }

    // Continuous at the join, and the shape function is exactly 1 there.
    assert.equal(elevToHeightKm(1, true), RELIEF_KM_PER_UNIT);
    assert.ok(Math.abs(elevToHeightKm(1 - 1e-9, true) - RELIEF_KM_PER_UNIT) < 1e-6);

    // Above the domain the shape is its own argument, so the branch is the
    // curve's scale times the parameter and introduces no new constant.
    for (const t of [1.0, 1.25, 1.5, 1.4280]) {
        assert.ok(Math.abs(elevToHeightKm(t, true) - RELIEF_KM_PER_UNIT * t) < 1e-9);
    }

    // Below the join nothing moved: the published shape values are unchanged.
    assert.ok(Math.abs(elevToHeightKm(0.5, true) - 1.125) < 1e-9);
    assert.ok(Math.abs(elevToHeightKm(0.75, true) - 3.7968750) < 1e-6);

    // The clamp's own failure mode, asserted so it cannot come back: the raw
    // polynomial is below the join value everywhere above it, and below sea
    // level past t = 1.25.
    const raw = t => RELIEF_KM_PER_UNIT * t * t * t * t * (5 - 4 * t);
    assert.ok(raw(1.1) < RELIEF_KM_PER_UNIT);
    assert.ok(raw(1.3) < 0);
});

test('surface classes use ocean connectivity, not the sign of the elevation', () => {
    const elev = new Float32Array([-1, -0.5, 0.5]);
    const isOpenOcean = new Uint8Array([1, 0, 0]);   // cell 1 is a closed below-sea-level floor
    const inland = new Uint8Array([0, 1, 0]);
    const sc = surfaceClasses(elev, isOpenOcean, inland);
    assert.equal(sc[0], 0, 'connected sea is ocean');
    assert.equal(sc[1], 2, 'flooded closed basin is inland water, not ocean');
    assert.equal(surfaceClasses(elev, isOpenOcean, new Uint8Array(3))[1], 1,
        'dry closed basin below sea level is land');
});

// ── Final-terrain drainage ───────────────────────────────────────────────

test('final catchments resolve to terminals and cover the basin floor', () => {
    const c = ctx();
    assert.ok(c.basins.selected.length > 0, 'fixture produced no preserved basins');
    for (let i = 0; i < c.basins.finalCatchments.length; i++) {
        const cat = c.basins.finalCatchments[i];
        assert.ok(cat.catchmentAreaKm2 > 0, `basin ${cat.id} has an empty catchment`);
        assert.ok(Number.isFinite(cat.catchmentToFloorRatio));
    }
});

test('THE PRODUCTS AGREE: summing drainage_terminal reproduces every catchment area', () => {
    const c = ctx();
    const audit = c.basins.drainageConsistency;
    assert.ok(audit, 'the export must publish the consistency check, not just pass it');
    assert.equal(audit.disagreeing, 0,
        `${audit.disagreeing} basins where the field and the table disagree`);
    assert.ok(audit.worstRelativeError < 1e-6);
    // The two totals a consumer would compute independently.
    assert.ok(Math.abs(audit.totalByFieldKm2 - audit.totalDeclaredKm2)
        < 1e-6 * audit.totalDeclaredKm2);
});

test('the endorheic fraction names its denominator and uses surface_class', () => {
    const c = ctx();
    const d = c.basins.drainageConsistency;
    assert.equal(d.landDenominator, 'surface_class == land (subaerial)');
    let byClass = 0, bySign = 0;
    for (let r = 0; r < c.mesh.numRegions; r++) {
        if (c.basins.r_surfaceClass[r] === 1) byClass += c.cellArea[r];
        if (c.r_elevation[r] > 0) bySign += c.cellArea[r];
    }
    assert.ok(Math.abs(d.landAreaKm2 - byClass) < 1e-6 * byClass);
    assert.ok(d.landAreaKm2 >= bySign, 'subaerial land must not be the smaller set');
    assert.ok(Math.abs(d.fractionOfLand - d.totalDeclaredKm2 / byClass) < 1e-9);
    // The land_mask denominator is the trap; it must give a larger number.
    if (byClass > bySign) {
        assert.ok(d.totalDeclaredKm2 / bySign > d.fractionOfLand,
            'land_mask inflates the fraction — that is why the denominator is named');
    }
});

test('catchments root at the FINAL low point, not the catalogue sink', () => {
    const c = ctx();
    const term = c.basins.r_finalTerminal;
    let moved = 0;
    for (let i = 0; i < c.basins.finalCatchments.length; i++) {
        const cat = c.basins.finalCatchments[i];
        const fin = c.basins.finalMeasurements[i];
        assert.equal(cat.sink, fin.sink, 'catchment root must be the remeasured sink');
        assert.equal(term[cat.sink], cat.sink, 'every root must be its own terminal');
        if (cat.sink !== cat.catalogueSink) {
            moved++;
            // The whole point: the catalogue sink is no longer the low point,
            // so rooting there would put the terminal above its own catchment.
            assert.ok(c.r_elevation[cat.sink] <= c.r_elevation[cat.catalogueSink]);
        }
    }
    // Not asserting moved > 0 — a small fixture may erode nothing — but if any
    // moved, the catalogue sink must not be a root.
    for (const cat of c.basins.finalCatchments) {
        if (cat.sink !== cat.catalogueSink) {
            assert.notEqual(term[cat.catalogueSink], cat.catalogueSink);
        }
    }
    assert.ok(moved >= 0);
});

test('REGRESSION: rooting at the catalogue sink after the low point moves', () => {
    // Reproduce the reported contradiction directly. Erosion deepens a cell
    // that is not the catalogued sink, so the basin's low point moves; the
    // per-basin record then names the new cell while the routing still roots
    // at the old one.
    const w = craterWorld();
    const basins = detectBasins(w.mesh, w.elev, w.isOcean, {});
    const crater = basins.reduce((a, b) => (b.volumeKm3 > a.volumeKm3 ? b : a));
    const area = uniformCellArea(w.mesh.numRegions);

    const eroded = Float32Array.from(w.elev);
    const moved = Array.from(crater.members).find(r => r !== crater.sink);
    assert.ok(moved !== undefined);
    eroded[moved] = w.elev[crater.sink] - 0.05;            // new deepest cell

    const isOcean = new Uint8Array(w.mesh.numRegions);
    for (let r = 0; r < eroded.length; r++) if (eroded[r] <= 0) isOcean[r] = 1;
    const remeasured = remeasureBasins(w.mesh, eroded, isOcean, [crater], area);
    assert.equal(remeasured[0].sink, moved, 'remeasure should follow the low point');

    // The bug: route at the catalogue sink, publish the remeasured one.
    const stale = computeFinalPreservedCatchments(w.mesh, eroded, isOcean, [crater], area);
    const staleAudit = auditCatchmentConsistency(stale[0].terminal,
        [{ sink: remeasured[0].sink, catchmentAreaKm2: stale[0].catchmentAreaKm2 }], area);
    assert.equal(staleAudit.disagreeing, 1, 'this is the shipped contradiction');
    assert.ok(staleAudit.totalByFieldKm2 < staleAudit.totalDeclaredKm2,
        'the field accumulates nothing at the published sink');

    // The fix: route at the sink the record names.
    const fixed = computeFinalPreservedCatchments(w.mesh, eroded, isOcean, [crater], area,
        { sinks: [remeasured[0].sink] });
    assert.equal(fixed[0].sink, remeasured[0].sink);
    assert.equal(fixed[0].catalogueSink, crater.sink);
    const audit = auditCatchmentConsistency(fixed[0].terminal, fixed, area);
    assert.equal(audit.disagreeing, 0);
    assert.ok(audit.agrees);
});

test('a catchment root outside its own basin is rejected', () => {
    const w = craterWorld();
    const basins = detectBasins(w.mesh, w.elev, w.isOcean, {});
    const crater = basins.reduce((a, b) => (b.volumeKm3 > a.volumeKm3 ? b : a));
    const outside = Array.from({ length: w.mesh.numRegions }, (_, r) => r)
        .find(r => !Array.from(crater.members).includes(r));
    assert.throws(
        () => computeFinalPreservedCatchments(w.mesh, w.elev, w.isOcean, [crater], null,
            { sinks: [outside] }),
        /not one of its members/);
});

test('flow accumulation runs to the basin sink, not to the rim', () => {
    const w = craterWorld();
    const basins = detectBasins(w.mesh, w.elev, w.isOcean, {});
    const crater = basins.reduce((a, b) => (b.volumeKm3 > a.volumeKm3 ? b : a));
    const area = uniformCellArea(w.mesh.numRegions);

    const sinkMask = new Uint8Array(w.mesh.numRegions);
    sinkMask[crater.sink] = 1;
    const bySink = computeRivers(w.mesh, w.elev, w.isOcean, area, { terminalRoot: sinkMask });

    // The old wiring: every member cell terminates flow.
    const memberMask = new Uint8Array(w.mesh.numRegions);
    for (const r of crater.members) memberMask[r] = 1;
    const byMember = computeRivers(w.mesh, w.elev, w.isOcean, area, { terminalRoot: memberMask });

    assert.ok(bySink.flowAccumKm2[crater.sink] > byMember.flowAccumKm2[crater.sink],
        'terminating at the member mask starves the sink of its own inflow');
    assert.ok(bySink.flowAccumKm2[crater.sink] >= crater.members.length * area[0] * 0.9,
        'the sink should carry at least the basin floor it drains');
});

test('flow accumulation conserves area over the whole tree', () => {
    const w = craterWorld();
    const area = uniformCellArea(w.mesh.numRegions);
    const { flowAccumKm2, drainTo } = computeRivers(w.mesh, w.elev, w.isOcean, area);
    let outflow = 0;
    for (let r = 0; r < w.mesh.numRegions; r++) if (drainTo[r] < 0) outflow += flowAccumKm2[r];
    let total = 0;
    for (let r = 0; r < w.mesh.numRegions; r++) total += area[r];
    assert.ok(Math.abs(outflow - total) < 1e-3 * total,
        `roots carry ${outflow} of ${total} — the tree leaks or double-counts`);
});

test('compareCatchments is a correct Jaccard index', () => {
    assert.equal(compareCatchments([1, 2, 3], [1, 2, 3], 10).jaccard, 1);
    assert.equal(compareCatchments([1, 2], [3, 4], 10).jaccard, 0);
    const partial = compareCatchments([1, 2, 3], [2, 3, 4], 10);
    assert.equal(partial.sharedCount, 2);
    assert.equal(partial.jaccard, 2 / 4);
    assert.equal(partial.addedCount, 1);
    assert.equal(partial.removedCount, 1);
    assert.equal(compareCatchments([], [], 10).jaccard, 1);
});

test('drainage terminals are open ocean or preserved sinks, never a mid-slope cell', () => {
    const w = craterWorld();
    const basins = detectBasins(w.mesh, w.elev, w.isOcean, {});
    const crater = basins.reduce((a, b) => (b.volumeKm3 > a.volumeKm3 ? b : a));
    const roots = new Uint8Array(w.mesh.numRegions);
    roots[crater.sink] = 1;
    const { terminal } = computeFinalTerrainDrainage(w.mesh, w.elev, w.isOcean, roots);
    const { isOpenOcean } = findOpenOcean(w.mesh, w.isOcean);
    for (let r = 0; r < w.mesh.numRegions; r++) {
        const t = terminal[r];
        if (t < 0) continue;
        assert.ok(isOpenOcean[t] || roots[t], `cell ${r} drains to ${t}, which is not a terminal`);
    }
});

// ── End-to-end ───────────────────────────────────────────────────────────

test('end to end: preservation beats vanilla drainage on the same depressions', () => {
    const c = ctx();
    const off = quiet(() => runGeneratePipeline({ ...PARAMS, preserveBasins: false }));
    const sel = c.basins.selected;
    const oceanOf = (e) => { const m = new Uint8Array(e.length); for (let r = 0; r < e.length; r++) if (e[r] <= 0) m[r] = 1; return m; };
    const vanilla = remeasureBasins(off.mesh, off.r_elevation, oceanOf(off.r_elevation), sel, c.cellArea);

    const weight = (arr, pick) => {
        let w = 0, v = 0;
        for (let i = 0; i < sel.length; i++) { w += sel[i].areaKm2; v += sel[i].areaKm2 * (pick(arr[i]) ?? 0); }
        return v / w;
    };
    const kept = weight(c.basins.finalMeasurements, m => m.retainedFraction);
    const lost = weight(vanilla, m => m.retainedFraction);

    assert.ok(kept > 0.5, `preserved retention ${(100 * kept).toFixed(1)}% should clear 50%`);
    assert.ok(lost < 0.15, `vanilla retention ${(100 * lost).toFixed(1)}% should be near zero`);
    assert.ok(kept > lost * 5, 'preservation must be dramatically better than vanilla');
});

test('preserveBasins:false leaves no basin state at all', () => {
    const off = quiet(() => runGeneratePipeline({ ...PARAMS, preserveBasins: false }));
    assert.equal(off.basins, null);
});

test('below-sea-level closed basins survive fixupTopology when resolvable', () => {
    const c = ctx();
    assert.ok(c.basins.topologyFixup, 'fixupTopology reported nothing');
    assert.ok(c.basins.topologyFixup.keptComponents > 0,
        'every below-sea-level closed basin was filled — fixupTopology is still destroying them');
});

// ── Drainage hypothesis (preserve / carve lists) ─────────────────────────

test('parseBasinList reads text, JSON, comments and retain fractions', () => {
    const txt = parseBasinList(`
        # a comment
        orogen-basin-aaa
        orogen-basin-bbb   0.35
        orogen-basin-ccc,0     # trailing comment
    `);
    assert.deepEqual(txt, [
        { id: 'orogen-basin-aaa', retain: 1 },
        { id: 'orogen-basin-bbb', retain: 0.35 },
        { id: 'orogen-basin-ccc', retain: 0 },
    ]);
    assert.deepEqual(parseBasinList('["a","b"]'), [{ id: 'a', retain: 1 }, { id: 'b', retain: 1 }]);
    assert.deepEqual(parseBasinList('[{"id":"a","retain":0.5}]'), [{ id: 'a', retain: 0.5 }]);
    assert.throws(() => parseBasinList('a 1.5'), /must be in \[0,1\]/);
    assert.throws(() => parseBasinList('a -0.1'), /must be in \[0,1\]/);
});

test('THE LOOP INVARIANT: basin ids do not depend on the carve decision', () => {
    // The whole iterative water-balance loop rests on this. Detection runs on
    // the pre-conditioning surface, which is produced before any preserve or
    // carve decision, so a verdict list computed against one export stays valid
    // against the next.
    const runs = [
        quiet(() => runGeneratePipeline({ ...PARAMS })),
        quiet(() => runGeneratePipeline({ ...PARAMS, basinMinCells: 60 })),
        quiet(() => runGeneratePipeline({ ...PARAMS, basinMinCells: 400 })),
    ];
    const ids = runs.map(c => c.basins.catalogue.map(b => b.id));
    assert.ok(runs[0].basins.selected.length !== runs[2].basins.selected.length,
        'the fixture must actually vary the preserved set for this to mean anything');
    assert.deepEqual(ids[0], ids[1]);
    assert.deepEqual(ids[0], ids[2]);
    assert.equal(hashJson(ids[0]), hashJson(ids[2]), 'catalogue hash must be identical');
});

test('a preserve list replaces threshold selection exactly', () => {
    const base = quiet(() => runGeneratePipeline({ ...PARAMS }));
    const wanted = base.basins.selected.slice(0, 5).map(b => ({ id: b.id, retain: 1 }));
    const c = quiet(() => runGeneratePipeline({ ...PARAMS, preserveBasinList: wanted }));
    assert.equal(c.basins.selected.length, 5);
    assert.deepEqual(c.basins.selected.map(b => b.id).sort(), wanted.map(w => w.id).sort());
    assert.equal(c.basins.selectionSource, 'preserve-list');
});

test('a carve list subtracts from threshold selection', () => {
    const base = quiet(() => runGeneratePipeline({ ...PARAMS }));
    const drop = base.basins.selected.slice(0, 3).map(b => ({ id: b.id, retain: 0 }));
    const c = quiet(() => runGeneratePipeline({ ...PARAMS, carveBasinList: drop }));
    assert.equal(c.basins.selected.length, base.basins.selected.length - 3);
    const kept = new Set(c.basins.selected.map(b => b.id));
    for (const d of drop) assert.ok(!kept.has(d.id), `${d.id} should have been carved`);
});

test('an unknown id in a list warns rather than failing silently', () => {
    const c = quiet(() => runGeneratePipeline({
        ...PARAMS, preserveBasinList: [{ id: 'orogen-basin-nope', retain: 1 }],
    }));
    assert.ok(c.basins.warnings.some(w => /not in this catalogue/.test(w)));
});

test('retain 0 in a preserve list is a carve verdict, not a full-depth notch', () => {
    // Downstream sends a COMPLETE verdict — every basin listed, carved ones at
    // retain 0 rather than omitted. Selecting those would make each one a
    // terminal drainage root with an unprotected rim: water routed inward to a
    // sink the carve is free to breach, which is not a state terrain can mean.
    const base = quiet(() => runGeneratePipeline({ ...PARAMS }));
    const all = base.basins.selected.map(b => b.id);
    assert.ok(all.length >= 4);
    const keep = all.slice(0, 2), drop = all.slice(2);

    const c = quiet(() => runGeneratePipeline({
        ...PARAMS,
        preserveBasinList: [
            ...keep.map(id => ({ id, retain: 1 })),
            ...drop.map(id => ({ id, retain: 0 })),
        ],
    }));
    const got = c.basins.selected.map(b => b.id);
    assert.deepEqual(got.slice().sort(), keep.slice().sort(),
        'only the retain > 0 entries may be preserved');
    assert.equal(c.basins.carvedByList, drop.length);
    for (const b of c.basins.selected) assert.ok(b.retain > 0);
    // And none of the carved ones is a drainage terminal.
    const sinks = new Set(c.basins.finalCatchments.map(x => x.sink));
    for (const id of drop) {
        const b = base.basins.selected.find(x => x.id === id);
        assert.ok(!sinks.has(b.sink), `carved basin ${id} is still a terminal root`);
    }
});

test('partial retain actually incises the rim, and does so monotonically', () => {
    const base = quiet(() => runGeneratePipeline({ ...PARAMS }));
    const pick = base.basins.selected.slice(0, 12).map(b => b.id);
    const depthAt = (retain) => {
        const c = quiet(() => runGeneratePipeline({
            ...PARAMS, preserveBasinList: pick.map(id => ({ id, retain })),
        }));
        const v = c.basins.finalMeasurements
            .map(m => m.retainedFraction).filter(x => x !== null).sort((a, b) => a - b);
        return v[Math.floor(v.length / 2)];
    };
    const full = depthAt(1), half = depthAt(0.45);
    assert.ok(half < full - 0.05,
        `retain 0.45 kept ${half.toFixed(2)} vs ${full.toFixed(2)} at retain 1 — `
        + 'partial incision is inert, which is worse than not offering it');
});

test('retain 1 is bit-identical to having no allowance at all', () => {
    // A knob must not move the baseline. The first version stored an absolute
    // pre-erosion floor, which shifted mean land elevation on every planet even
    // with no retain fraction anywhere.
    const w = craterWorld();
    const basins = detectBasins(w.mesh, w.elev, w.isOcean, {});
    const crater = basins.reduce((a, b) => (b.volumeKm3 > a.volumeKm3 ? b : a));
    const prot = buildBasinProtection(w.mesh, [crater], {});
    for (let r = 0; r < w.mesh.numRegions; r++) {
        assert.equal(prot.carveAllowance[r], 0, 'no retain given, so no allowance');
    }
    const withAllow = Float32Array.from(w.elev);
    const without = Float32Array.from(w.elev);
    const stripped = { ...prot, carveAllowance: null };
    priorityFloodCarve(w.mesh, withAllow, w.isOcean, 1, prot);
    priorityFloodCarve(w.mesh, without, w.isOcean, 1, stripped);
    assert.deepEqual(Array.from(withAllow), Array.from(without));
});

test('incision only touches basins with retain below 1', () => {
    const w = craterWorld();
    const basins = detectBasins(w.mesh, w.elev, w.isOcean, {});
    const crater = basins.reduce((a, b) => (b.volumeKm3 > a.volumeKm3 ? b : a));
    const untouched = Float32Array.from(w.elev);
    const r1 = inciseOutlets(w.mesh, untouched, [{ ...crater, retain: 1 }]);
    assert.equal(r1.basinsCut, 0);
    assert.deepEqual(Array.from(untouched), Array.from(w.elev));

    const cut = Float32Array.from(w.elev);
    const r2 = inciseOutlets(w.mesh, cut, [{ ...crater, retain: 0.5 }]);
    assert.equal(r2.basinsCut, 1);
    assert.ok(r2.cellsLowered > 0);
    for (let r = 0; r < cut.length; r++) {
        assert.ok(cut[r] <= w.elev[r] + 1e-9, `incision raised cell ${r}`);
    }
});

// ── Catchment routing ────────────────────────────────────────────────────

test('filled routing resolves every land cell; raw routing does not', () => {
    const c = ctx();
    const isOcean = new Uint8Array(c.mesh.numRegions);
    for (let r = 0; r < c.mesh.numRegions; r++) if (c.r_elevation[r] <= 0) isOcean[r] = 1;
    const roots = new Uint8Array(c.mesh.numRegions);
    for (const b of c.basins.selected) roots[b.sink] = 1;

    const unresolved = (raw) => {
        const { terminal } = computeFinalTerrainDrainage(c.mesh, c.r_elevation, isOcean, roots, { raw });
        let land = 0, bad = 0;
        for (let r = 0; r < c.mesh.numRegions; r++) {
            if (c.r_elevation[r] <= 0) continue;
            land++; if (terminal[r] < 0) bad++;
        }
        return bad / land;
    };
    assert.equal(unresolved(false), 0,
        'filled routing must resolve every land cell, or catchments undercount');
    assert.ok(unresolved(true) > 0.1,
        'the raw measure should strand a large share — that was the bug being fixed');
});

// ── Determinism / reproducibility ────────────────────────────────────────

test('the same seed produces identical basins and identical hashes', () => {
    const a = quiet(() => runGeneratePipeline({ ...PARAMS, preserveBasins: true }));
    const b = quiet(() => runGeneratePipeline({ ...PARAMS, preserveBasins: true }));
    assert.deepEqual(a.basins.selected.map(x => x.id), b.basins.selected.map(x => x.id));
    assert.equal(hashTypedArray(a.r_elevation), hashTypedArray(b.r_elevation));
    assert.equal(hashTypedArray(a.basins.r_basinIndex), hashTypedArray(b.basins.r_basinIndex));
    assert.equal(hashJson(a.basins.selected.map(x => x.id)), hashJson(b.basins.selected.map(x => x.id)));
});

// ── Gravity and the catalogue's kilometres ───────────────────────────────

test('basin ...Km fields carry the 1/g relief scaling, and reliefScale 1 is the baseline', () => {
    // These are published beside `elevation_km` in the same manifest under the
    // same units. They used to skip the 1/g scaling entirely, so on a 1.31-g
    // planet a basin depth and the surface height above it were 31% apart —
    // two different kilometres, and a water balance reads both.
    const w = craterWorld();
    const cellArea = uniformCellArea(w.mesh.numRegions);
    const at = (reliefScale) => detectBasins(w.mesh, w.elev, w.isOcean, { cellArea, reliefScale });

    const base = at(undefined);
    const unit = at(1);
    assert.deepEqual(unit.map(b => b.depthKm), base.map(b => b.depthKm),
        'reliefScale 1 must be bit-identical to omitting it — a knob may not move the baseline');

    const half = at(0.5);
    assert.equal(half.length, base.length, 'gravity must not change what depressions exist');
    let scaledPairs = 0;
    for (let i = 0; i < base.length; i++) {
        assert.equal(half[i].id, base[i].id, 'basin ids must be gravity-invariant');
        assert.equal(half[i].depth, base[i].depth, 'model-unit depth must not scale');
        // Only wholly-above-sea-level basins scale by exactly the factor:
        // negative heights do not scale, matching elevation_km's own rule.
        if (base[i].sinkElevationKm > 0) {
            assert.ok(Math.abs(half[i].depthKm / base[i].depthKm - 0.5) < 1e-6,
                `depthKm should halve at reliefScale 0.5: ${base[i].depthKm} -> ${half[i].depthKm}`);
            assert.ok(Math.abs(half[i].volumeKm3 / base[i].volumeKm3 - 0.5) < 1e-6,
                'volumeKm3 should halve too');
            scaledPairs++;
        }
    }
    assert.ok(scaledPairs > 0, 'no above-sea-level basin to check the scaling against');
});

test('THE INVARIANT: gravity moves no part of the preserved set', () => {
    // selectBasins runs during generation, inciseOutlets notches the rims of what
    // it returns and buildBasinProtection constrains erosion over them, so the
    // preserved set reaches r_elevation. The generator's stated invariant is that
    // gravity leaves the model terrain alone and enters at export -- erosion
    // constants are calibrated in model units, and pushing a scaled parameter
    // through the nonlinear height curve damps the effect instead of delivering
    // it. So neither WHICH basins are preserved nor the ORDER they are
    // enumerated in may move with g. The order is what basin_index depends on;
    // the membership is what the carve list is drawn from and what the terrain
    // is conditioned around. tools/test-integration.mjs holds the same invariant
    // at the level of the elevation hash; this one says which mechanism it is.
    //
    // Both keys are measured at reference gravity for exactly this reason:
    // orderingVolumeKm3 for the sort, selectionDepthKm for the depth floor. The
    // area and cell floors are horizontal and never had the problem.
    const base = ctx().basins.selected;
    assert.ok(base.length > 1, 'need at least two basins for an ordering to exist');
    const baseIds = base.map(b => b.id);

    for (const gravityMS2 of [0.5, 2, 4].map(f => f * EARTH.gravityMS2)) {
        const run = quiet(() => runGeneratePipeline({
            ...PARAMS, preserveBasins: true, planet: { gravityMS2 },
        })).basins.selected;
        assert.deepEqual(run.map(b => b.id), baseIds,
            `the preserved set and its order must not move at g=${gravityMS2}`);
    }

    // A CONTROL, so the equality above cannot pass vacuously. Comparing the floor
    // against the scaled `depthKm` -- the physically-tempting reading, and the one
    // this test exists to forbid -- would drop every basin whose depth on a
    // 2 g planet falls under the 50 m floor. There must be some, or this fixture
    // never exercises the difference and the assertions above are about nothing.
    const heavy = quiet(() => runGeneratePipeline({
        ...PARAMS, preserveBasins: true, planet: { gravityMS2: 2 * EARTH.gravityMS2 },
    })).basins.selected;
    const shallowOnHeavy = heavy.filter(b => b.depthKm < BASIN_MIN_DEPTH_KM);
    assert.ok(shallowOnHeavy.length > 0,
        'no preserved basin is physically shallower than the floor at 2 g, so this '
        + 'fixture cannot tell the two currencies apart');
    for (const b of shallowOnHeavy) {
        assert.ok(b.selectionDepthKm >= BASIN_MIN_DEPTH_KM,
            'every preserved basin must clear the floor in the currency it is compared in');
    }
});

test('sha256 matches known vectors', () => {
    assert.equal(sha256(new TextEncoder().encode('')),
        'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855');
    assert.equal(sha256(new TextEncoder().encode('abc')),
        'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad');
});

test('hashJson is insensitive to key order but not to values', () => {
    assert.equal(hashJson({ a: 1, b: 2 }), hashJson({ b: 2, a: 1 }));
    assert.notEqual(hashJson({ a: 1 }), hashJson({ a: 2 }));
});
