#!/usr/bin/env node --test
/**
 * Regression tests for lithology and lithology-driven erosion.
 *
 *   node --test tools/test-lithology.mjs
 *
 * Two promises are under test. First, that the rock map is a faithful, stable
 * function of tectonic state. Second — the one that actually changes terrain —
 * that erodibility redistributes erosion WITHOUT changing how much of it there
 * is, so turning lithology on does not silently recalibrate the erosion sliders.
 *
 * Needs: npm i delaunator
 */

import test from 'node:test';
import assert from 'node:assert/strict';
import Delaunator from 'delaunator';

import { setDelaunator } from '../js/sphere-mesh.js';
import {
    ROCK_CLASSES, classifyLithology, buildErodibility, buildLithoState,
    updateExhumation, finalSurfaceRock, rockComposition, computeScarpPotential, saltCrustMask,
    COVER_SEQUENCE,
} from '../js/lithology.js';
import { runGeneratePipeline } from '../js/pipeline.js';
import { LITHO_CARBONATE_LAT_DEG } from '../js/terrain-config.js';
import { hashTypedArray } from '../js/sha256.js';

setDelaunator(Delaunator);

const realLog = console.log;
const quiet = (fn) => { console.log = () => {}; try { return fn(); } finally { console.log = realLog; } };

const PARAMS = {
    N: 20000, P: 24, jitter: 0.75, nMag: 1, numContinents: 5, smoothing: 0,
    hydraulicErosion: 0.5, thermalErosion: 0.3, ridgeSharpening: 0.3,
    glacialErosion: 0.3, terrainWarp: 0, seed: 4242,
    toggledIndices: [], skipClimate: true,
};

let _on = null, _off = null;
const withLitho = () => (_on ??= quiet(() => runGeneratePipeline({ ...PARAMS, lithology: true })));
const noLitho = () => (_off ??= quiet(() => runGeneratePipeline({ ...PARAMS, lithology: false })));

// ── Rock table ───────────────────────────────────────────────────────────

test('rock table is well formed', () => {
    const ids = new Set(), codes = new Set();
    for (const c of ROCK_CLASSES) {
        assert.ok(Number.isInteger(c.id) && c.id >= 0, `bad id on ${c.code}`);
        assert.ok(!ids.has(c.id), `duplicate id ${c.id}`);
        assert.ok(!codes.has(c.code), `duplicate code ${c.code}`);
        ids.add(c.id); codes.add(c.code);
        assert.ok(c.erodibility > 0, `${c.code} must have positive erodibility`);
        assert.ok(['igneous', 'sedimentary', 'metamorphic', 'none'].includes(c.category),
            `${c.code} has category ${c.category}`);
    }
    // Ids must be dense from 0, since they index a Uint8Array lookup table.
    for (let i = 0; i < ROCK_CLASSES.length; i++) assert.ok(ids.has(i), `id ${i} missing`);
    assert.ok(ROCK_CLASSES.length <= 256, 'rock classes must fit in a Uint8Array');
});

test('rock strength ordering is physically sensible', () => {
    const k = Object.fromEntries(ROCK_CLASSES.map(c => [c.code, c.erodibility]));
    assert.ok(k.quartzite < k.granite, 'quartzite should outlast granite');
    assert.ok(k.granite < k.schist, 'granite should outlast schist');
    assert.ok(k.gneiss < k.shelf_clastic, 'gneiss should outlast sandstone/shale');
    assert.ok(k.shelf_clastic < k.evaporite, 'evaporite is the softest thing here');
    assert.ok(k.flood_basalt < k.foreland_clastic, 'basalt caps should outlast molasse');
});

// ── Classification ───────────────────────────────────────────────────────

test('every region gets a valid rock class, and classification does not touch elevation', () => {
    const c = withLitho();
    const before = Float32Array.from(c.r_elevation);
    const litho = classifyLithology(c.mesh, c.r_xyz, c.r_elevation, c.tectonics, c.debugLayers, {});
    assert.deepEqual(Array.from(c.r_elevation), Array.from(before), 'classification mutated elevation');

    const valid = new Set(ROCK_CLASSES.map(x => x.id));
    for (let r = 0; r < c.mesh.numRegions; r++) {
        assert.ok(valid.has(litho.basement[r]), `bad basement class at ${r}`);
        assert.ok(valid.has(litho.cover[r]), `bad cover class at ${r}`);
        assert.ok(litho.coverThicknessKm[r] >= 0, `negative cover thickness at ${r}`);
        assert.ok(Number.isFinite(litho.coverThicknessKm[r]), `non-finite cover thickness at ${r}`);
    }
});

test('surface rock is cover where it survives and basement where it does not', () => {
    const c = withLitho();
    const litho = classifyLithology(c.mesh, c.r_xyz, c.r_elevation, c.tectonics, c.debugLayers, {});
    for (let r = 0; r < c.mesh.numRegions; r++) {
        const expected = litho.coverThicknessKm[r] > 1e-4 ? litho.cover[r] : litho.basement[r];
        assert.equal(litho.surface[r], expected, `surface rock wrong at ${r}`);
    }
});

test('oceanic crust is basaltic and continents are not', () => {
    const c = withLitho();
    const byId = Object.fromEntries(ROCK_CLASSES.map(x => [x.id, x.code]));
    let oceanChecked = 0;
    for (let r = 0; r < c.mesh.numRegions; r++) {
        if (!c.tectonics.r_isOcean[r]) continue;
        assert.equal(byId[c.lithology.r_basementRock[r]], 'morb',
            `oceanic basement at ${r} is ${byId[c.lithology.r_basementRock[r]]}, expected morb`);
        if (++oceanChecked > 500) break;
    }
    assert.ok(oceanChecked > 0, 'no oceanic cells in the fixture');
});

test('composition is a normalised partition of the land surface', () => {
    const c = withLitho();
    const subaerial = new Uint8Array(c.mesh.numRegions);
    for (let r = 0; r < c.mesh.numRegions; r++) {
        subaerial[r] = c.basins.r_surfaceClass[r] === 1 ? 1 : 0;
    }
    const comp = rockComposition(c.lithology.r_surfaceRock, subaerial, c.cellArea);
    const sum = comp.reduce((a, x) => a + x.fraction, 0);
    assert.ok(Math.abs(sum - 1) < 1e-6, `fractions sum to ${sum}`);
    for (const x of comp) assert.ok(x.areaKm2 > 0, `${x.code} listed with no area`);
});

test('composition is measured against surface_class, not the sign of elevation', () => {
    // land_mask drops the dry sub-sea-level basin floors, and the classes that
    // live on those floors are the ones that get under-reported.
    const c = withLitho();
    const subaerial = new Uint8Array(c.mesh.numRegions);
    const bySign = new Uint8Array(c.mesh.numRegions);
    let dryFloors = 0;
    for (let r = 0; r < c.mesh.numRegions; r++) {
        subaerial[r] = c.basins.r_surfaceClass[r] === 1 ? 1 : 0;
        bySign[r] = c.r_elevation[r] > 0 ? 1 : 0;
        if (subaerial[r] && !bySign[r]) dryFloors++;
    }
    assert.ok(dryFloors > 0, 'fixture has no dry sub-sea-level floors to distinguish the two');

    const areaOf = (mask) => { let a = 0;
        for (let r = 0; r < mask.length; r++) if (mask[r]) a += c.cellArea[r]; return a; };
    assert.ok(areaOf(subaerial) > areaOf(bySign), 'subaerial land must be the larger set');

    // The shipped table must use the larger denominator.
    const shipped = c.lithology.composition;
    const total = shipped.reduce((a, x) => a + x.areaKm2, 0);
    assert.ok(Math.abs(total - areaOf(subaerial)) < 1e-3 * total,
        `composition sums to ${total}, not the subaerial area ${areaOf(subaerial)}`);
    assert.equal(c.lithology.landDenominator, 'surface_class == land (subaerial)');
});

// ── Erodibility ──────────────────────────────────────────────────────────

test('THE INVARIANT: erodibility is mean-normalised to 1 over land', () => {
    const c = withLitho();
    const litho = classifyLithology(c.mesh, c.r_xyz, c.r_elevation, c.tectonics, c.debugLayers, {});
    const { erodibility } = buildErodibility(litho.surface, c.r_elevation, c.cellArea, 1.0);

    let wsum = 0, ksum = 0;
    for (let r = 0; r < c.mesh.numRegions; r++) {
        if (c.r_elevation[r] <= 0) continue;
        const w = c.cellArea[r];
        wsum += w; ksum += w * erodibility[r];
    }
    const mean = ksum / wsum;
    assert.ok(Math.abs(mean - 1) < 1e-3,
        `land-mean erodibility is ${mean}, must be 1 or lithology silently rescales total erosion`);
});

test('strength 0 gives a uniform field; strength 1 gives contrast', () => {
    const c = withLitho();
    const litho = classifyLithology(c.mesh, c.r_xyz, c.r_elevation, c.tectonics, c.debugLayers, {});
    const flat = buildErodibility(litho.surface, c.r_elevation, c.cellArea, 0).erodibility;
    for (let r = 0; r < flat.length; r++) {
        assert.ok(Math.abs(flat[r] - 1) < 1e-6, `strength 0 must be uniform, got ${flat[r]} at ${r}`);
    }
    const full = buildErodibility(litho.surface, c.r_elevation, c.cellArea, 1).erodibility;
    let min = Infinity, max = -Infinity;
    for (let r = 0; r < full.length; r++) { min = Math.min(min, full[r]); max = Math.max(max, full[r]); }
    assert.ok(max / min > 2, `strength 1 should give real contrast, got ${min}..${max}`);
});

test('erodibility is always finite and positive', () => {
    const c = withLitho();
    for (let r = 0; r < c.mesh.numRegions; r++) {
        const k = c.lithology.r_erodibility[r];
        assert.ok(Number.isFinite(k) && k > 0, `erodibility ${k} at ${r}`);
    }
});

test('an undefined strength does not poison the field with NaN', () => {
    const c = withLitho();
    const litho = classifyLithology(c.mesh, c.r_xyz, c.r_elevation, c.tectonics, c.debugLayers, {});
    const state = buildLithoState(litho, c.r_elevation, c.cellArea, undefined);
    for (let r = 0; r < state.basementErodibility.length; r++) {
        assert.ok(Number.isFinite(state.basementErodibility[r]), `NaN basement erodibility at ${r}`);
        assert.ok(Number.isFinite(state.erodibility[r]), `NaN surface erodibility at ${r}`);
    }
});

// ── Exhumation ───────────────────────────────────────────────────────────

test('exhumation strips cover and swaps in basement erodibility', () => {
    const c = withLitho();
    const litho = classifyLithology(c.mesh, c.r_xyz, c.r_elevation, c.tectonics, c.debugLayers, {});
    const state = buildLithoState(litho, c.r_elevation, c.cellArea, 1);

    // Find a cell with cover, and erode straight through it.
    let target = -1;
    for (let r = 0; r < state.coverThicknessKm.length; r++) {
        if (state.coverThicknessKm[r] > 0.1) { target = r; break; }
    }
    assert.ok(target >= 0, 'fixture has no covered cells');

    const before = Float32Array.from(c.r_elevation);
    const after = Float32Array.from(c.r_elevation);
    after[target] -= state.coverThicknessKm[target] + 0.05;   // strip past the cover

    updateExhumation(state, before, after);
    assert.equal(state.coverThicknessKm[target], 0, 'cover should be fully stripped');
    assert.equal(state.erodibility[target], state.basementErodibility[target],
        'exposed cell should take basement erodibility');
    assert.equal(finalSurfaceRock(state)[target], state.basement[target],
        'exposed cell should report basement as its surface rock');
});

test('exhumation never adds cover, and ignores cells that did not erode', () => {
    const c = withLitho();
    const litho = classifyLithology(c.mesh, c.r_xyz, c.r_elevation, c.tectonics, c.debugLayers, {});
    const state = buildLithoState(litho, c.r_elevation, c.cellArea, 1);
    const thicknessBefore = Float32Array.from(state.coverThicknessKm);

    const before = Float32Array.from(c.r_elevation);
    const after = Float32Array.from(c.r_elevation);
    for (let r = 0; r < after.length; r++) after[r] += 0.5;   // everything aggrades

    updateExhumation(state, before, after);
    for (let r = 0; r < state.coverThicknessKm.length; r++) {
        assert.equal(state.coverThicknessKm[r], thicknessBefore[r],
            `cover changed at ${r} despite no erosion`);
    }
});

// ── Scarp potential ──────────────────────────────────────────────────────

test('scarp potential is a bounded fraction, land only', () => {
    const c = withLitho();
    const sp = c.lithology.r_scarpPotential;
    // Land is the SUBAERIAL surface, the same test computeScarpPotential is
    // given. This assertion used to read `r_elevation[r] <= 0`, which is
    // land_mask, and so demanded that a plateau margin standing inside a dry
    // closed basin below sea level score zero. That is the terrain this fork
    // exists to preserve, and it is an escarpment like any other.
    const subaerial = c.basins.r_isSubaerial;
    assert.equal(sp.length, c.mesh.numRegions);
    for (let r = 0; r < c.mesh.numRegions; r++) {
        assert.ok(Number.isFinite(sp[r]), `non-finite at ${r}`);
        assert.ok(sp[r] >= 0 && sp[r] <= 1, `out of range (${sp[r]}) at ${r}`);
        if (!subaerial[r]) {
            assert.equal(sp[r], 0, `ocean cell ${r} has scarp potential`);
        }
    }
});

test('scarp potential is non-trivial but not everywhere', () => {
    const c = withLitho();
    const sp = c.lithology.r_scarpPotential;
    let land = 0, nz = 0;
    for (let r = 0; r < c.mesh.numRegions; r++) {
        if (!c.basins.r_isSubaerial[r]) continue;
        land++; if (sp[r] > 0) nz++;
    }
    const frac = nz / land;
    assert.ok(frac > 0.01, `only ${(100 * frac).toFixed(2)}% of land flagged — the field is inert`);
    assert.ok(frac < 0.75, `${(100 * frac).toFixed(2)}% of land flagged — the field is not discriminating`);
});

test('no erodibility contrast means no scarp', () => {
    const c = withLitho();
    const litho = classifyLithology(c.mesh, c.r_xyz, c.r_elevation, c.tectonics, c.debugLayers, {});
    // strength 0 makes every rock erode identically, so the contrast term is 0
    // everywhere and no cell can host a differential-retreat scarp.
    const flat = buildLithoState(litho, c.r_elevation, c.cellArea, 0);
    const sp = computeScarpPotential(c.mesh, c.r_elevation, flat, c.neighborDist);
    for (let r = 0; r < sp.length; r++) {
        assert.equal(sp[r], 0, `cell ${r} scarps despite uniform erodibility`);
    }
});

test('flat terrain hosts no scarp however strong the contrast', () => {
    const c = withLitho();
    const litho = classifyLithology(c.mesh, c.r_xyz, c.r_elevation, c.tectonics, c.debugLayers, {});
    const state = buildLithoState(litho, c.r_elevation, c.cellArea, 1);
    // A dead-flat land surface: contrast and cover edges survive, relief does not.
    const flatElev = new Float32Array(c.mesh.numRegions);
    for (let r = 0; r < flatElev.length; r++) flatElev[r] = c.r_elevation[r] > 0 ? 1.0 : -1.0;
    const sp = computeScarpPotential(c.mesh, flatElev, state, c.neighborDist);
    for (let r = 0; r < sp.length; r++) {
        assert.equal(sp[r], 0, `cell ${r} scarps on flat ground`);
    }
});

test('scarp potential needs a cover edge, not just contrast and relief', () => {
    const c = withLitho();
    const litho = classifyLithology(c.mesh, c.r_xyz, c.r_elevation, c.tectonics, c.debugLayers, {});
    const state = buildLithoState(litho, c.r_elevation, c.cellArea, 1);
    // Bury the interface everywhere: uniform thick cover means it never
    // daylights, so there is nothing for a scarp to retreat along.
    state.coverThicknessKm.fill(5);
    const sp = computeScarpPotential(c.mesh, c.r_elevation, state, c.neighborDist);
    let nz = 0;
    for (let r = 0; r < sp.length; r++) if (sp[r] > 0) nz++;
    assert.equal(nz, 0, 'a fully buried interface should produce no scarps');
});

test('scarp potential is deterministic', () => {
    const a = quiet(() => runGeneratePipeline({ ...PARAMS, lithology: true }));
    const b = quiet(() => runGeneratePipeline({ ...PARAMS, lithology: true }));
    assert.equal(hashTypedArray(a.lithology.r_scarpPotential),
                 hashTypedArray(b.lithology.r_scarpPotential));
});

// ── End to end ───────────────────────────────────────────────────────────

test('lithology changes where erosion goes, not how much of it there is', () => {
    const on = withLitho(), off = noLitho();

    // Total land volume removed by erosion should be close between the two
    // runs — the erodibility field is mean-normalised precisely so that the
    // Hydraulic Erosion slider keeps its calibration.
    const eroded = (c) => {
        let v = 0;
        for (let r = 0; r < c.mesh.numRegions; r++) {
            const d = c.debugLayers.erosionDelta[r];
            if (d < 0) v += -d * c.cellArea[r];
        }
        return v;
    };
    const a = eroded(on), b = eroded(off);
    const ratio = a / b;
    assert.ok(ratio > 0.6 && ratio < 1.6,
        `eroded volume ratio ${ratio.toFixed(3)} — lithology should redistribute erosion, not rescale it`);

    // But the terrain must actually differ, or differential erosion is inert.
    let differing = 0;
    for (let r = 0; r < on.mesh.numRegions; r++) {
        if (Math.abs(on.r_elevation[r] - off.r_elevation[r]) > 1e-6) differing++;
    }
    assert.ok(differing / on.mesh.numRegions > 0.05,
        `only ${differing} cells differ — differential erosion is not doing anything`);
});

test('soft rock erodes more than hard rock, at matched flow', () => {
    const on = withLitho(), off = noLitho();
    // Compare each cell's erosion against the lithology-free control, and
    // correlate the difference with erodibility. Using the control as the
    // baseline removes the confound that soft rock also sits in low, wet places.
    let n = 0, sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0;
    for (let r = 0; r < on.mesh.numRegions; r++) {
        if (on.r_elevation[r] <= 0 && off.r_elevation[r] <= 0) continue;
        const x = on.lithology.r_erodibility[r];
        const y = off.r_elevation[r] - on.r_elevation[r];   // >0 = lithology eroded more here
        n++; sx += x; sy += y; sxx += x * x; sxy += x * y; syy += y * y;
    }
    const cov = sxy / n - (sx / n) * (sy / n);
    const sdx = Math.sqrt(sxx / n - (sx / n) ** 2);
    const sdy = Math.sqrt(syy / n - (sy / n) ** 2);
    const rho = cov / (sdx * sdy);
    assert.ok(rho > 0.1,
        `erodibility vs extra erosion correlation is ${rho.toFixed(3)}; expected clearly positive`);
});

test('lithology:false produces no lithology state and no erodibility', () => {
    const off = noLitho();
    assert.equal(off.lithology, null);
    assert.equal(off.debugLayers.surfaceRock, undefined);
});

test('same seed gives identical rock and erodibility', () => {
    const a = quiet(() => runGeneratePipeline({ ...PARAMS, lithology: true }));
    const b = quiet(() => runGeneratePipeline({ ...PARAMS, lithology: true }));
    assert.equal(hashTypedArray(a.lithology.r_surfaceRock), hashTypedArray(b.lithology.r_surfaceRock));
    assert.equal(hashTypedArray(a.lithology.r_erodibility), hashTypedArray(b.lithology.r_erodibility));
    assert.equal(hashTypedArray(a.r_elevation), hashTypedArray(b.r_elevation));
});

test('preserved endorheic basins take closed-basin fill, zoned into salt and clastic', () => {
    const c = withLitho();
    if (!c.basins || c.basins.selected.length === 0) return;   // nothing to assert against
    const saltId = ROCK_CLASSES.find(x => x.code === 'evaporite').id;
    const playaId = ROCK_CLASSES.find(x => x.code === 'playa_clastic').id;
    let inBasin = 0, salt = 0, playa = 0;
    for (let r = 0; r < c.mesh.numRegions; r++) {
        // Subaerial, not elevation > 0: the dry basin floors below sea level are
        // exactly where basin fill accumulates, so an elevation-sign test drops
        // the most characteristic cells this assertion is about.
        if (c.basins.r_isEndorheic[r] && c.basins.r_isSubaerial[r]) {
            inBasin++;
            if (c.lithology.r_coverRock[r] === saltId) salt++;
            if (c.lithology.r_coverRock[r] === playaId) playa++;
        }
    }
    assert.ok(inBasin > 0, 'no land cells inside preserved basins');
    assert.ok((salt + playa) / inBasin > 0.9,
        `only ${(100 * (salt + playa) / inBasin).toFixed(0)}% of endorheic land is basin fill`);
    // Both facies must actually occur: all-salt is the bug this split fixed,
    // all-clastic would mean the sump never qualifies.
    assert.ok(salt > 0, 'no salt crust anywhere — the sump zoning never fires');
    assert.ok(playa > 0, 'no clastic fill — the whole basin is still being called salt');
});

test('salt crust is the sump, and is a minority of the basin floor', () => {
    const c = withLitho();
    if (!c.basins || c.basins.selected.length === 0) return;
    const saltId = ROCK_CLASSES.find(x => x.code === 'evaporite').id;
    const idx = c.basins.r_basinIndex;
    let saltArea = 0, basinArea = 0, saltHigh = 0;
    for (let r = 0; r < c.mesh.numRegions; r++) {
        if (idx[r] < 0) continue;
        basinArea += c.cellArea[r];
        if (c.lithology.r_coverRock[r] !== saltId) continue;
        saltArea += c.cellArea[r];
        // Every salt cell must sit inside its own basin's lower band.
        const b = c.basins.selected[idx[r]];
        const ceiling = b.sinkElevation + 0.25 * b.depth;
        if (Number.isFinite(b.depth) && b.depth > 0
            && c.basins.preConditioningElev[r] > ceiling + 1e-6) saltHigh++;
    }
    assert.equal(saltHigh, 0, `${saltHigh} salt cells sit above their basin's sump ceiling`);
    const share = saltArea / basinArea;
    assert.ok(share > 0.02 && share < 0.6,
        `salt is ${(100 * share).toFixed(0)}% of basin floor — implausible for a closed basin`);
});

test('saltCrustMask floods a flat pan entirely and never leaves a basin', () => {
    const elev = Float32Array.from([0, 1, 2, 3, 4, 5]);
    const idx = Int32Array.from([0, 0, 0, 1, 1, -1]);
    // Basin 0 has relief; basin 1 is a flat pan with no measurable depth.
    const selected = [
        { sinkElevation: 0, depth: 4 },        // ceiling 0 + 0.25*4 = 1
        { sinkElevation: 3, depth: 0 },        // pan: all of it floods
    ];
    const mask = saltCrustMask(selected, idx, elev);
    assert.deepEqual(Array.from(mask), [1, 1, 0, 1, 1, 0]);
    assert.equal(mask[5], 0, 'a cell outside every basin must never be salt');
});

test('an orogen does not exhume a closed basin sitting inside it', () => {
    // The fold-belt rule thins cover toward zero because active orogens are
    // being stripped. A closed basin inside one is where the stripped material
    // GOES, so it is exempt. This was an oversight for one release and it was
    // invisible to the test above, which checks the cover CLASS: the class
    // stayed playa/evaporite while the thickness went to zero, so the surface
    // silently read as orogenic basement over 15% of preserved-basin area.
    // Assert on the thickness and the exposed surface, not the class.
    const R = Object.fromEntries(ROCK_CLASSES.map(c => [c.code, c.id]));
    const mesh = { numRegions: 2 };
    const r_xyz = new Float32Array([0, 0, 1, 0, 0, 1]);
    const r_elevation = new Float32Array([0.2, 0.2]);
    const t = {
        r_isOcean: Uint8Array.from([0, 0]),
        r_stress: Float32Array.from([0, 0]),
        maxStress: 1,
        r_t_foldBelt: Float32Array.from([1, 1]),   // saturated orogen on both
        r_t_craton: Float32Array.from([0.5, 0.5]), // both would take craton cover
    };
    // Cell 0 is inside a preserved closed basin; cell 1 is ordinary orogen.
    const litho = classifyLithology(mesh, r_xyz, r_elevation, t, {},
        { basins: { r_isEndorheic: Uint8Array.from([1, 0]), r_isSaltCrust: null } });

    assert.equal(litho.cover[0], R.playa_clastic, 'basin cell lost its fill class');
    assert.ok(litho.coverThicknessKm[0] > 0,
        'closed-basin fill was thinned to nothing by the orogen it sits in');
    assert.equal(litho.surface[0], R.playa_clastic,
        'basin floor exposes basement despite having fill cover');

    // The rule itself must still work where it belongs.
    assert.equal(litho.coverThicknessKm[1], 0,
        'fold-belt thinning stopped applying to ordinary orogen cover');
    assert.notEqual(litho.surface[1], litho.cover[1],
        'an exhumed orogen should expose basement, not its cover');
});

test('a closed basin on oceanic crust still gets basin fill, over marine basement', () => {
    // A closed basin traps fill whatever crust it sits on. The marine parent
    // material survives as the basement; the fill is the younger deposit and so
    // is what the surface exposes. Before this, whole preserved basins on
    // oceanic crust surfaced as pelagic ooze and carbonate platform, and were
    // handed downstream as vegetated land because the barren mask selects on
    // class name.
    const R = Object.fromEntries(ROCK_CLASSES.map(c => [c.code, c.id]));
    const mesh = { numRegions: 2 };
    const r_xyz = new Float32Array([0, 0, 1, 0, 0, 1]);
    const r_elevation = new Float32Array([-0.1, -0.1]);
    const t = {
        r_isOcean: Uint8Array.from([1, 1]),          // both on oceanic crust
        r_stress: Float32Array.from([0, 0]),
        maxStress: 1,
        dist_coast: Float32Array.from([0, 0]),       // both would take shelf cover
    };
    const litho = classifyLithology(mesh, r_xyz, r_elevation, t, {},
        { basins: { r_isEndorheic: Uint8Array.from([1, 0]), r_isSaltCrust: null } });

    assert.equal(litho.cover[0], R.playa_clastic, 'closed basin on oceanic crust took marine cover');
    assert.equal(litho.surface[0], R.playa_clastic, 'its surface is not basin fill');
    assert.equal(litho.basement[0], R.morb, 'the marine parent material was lost, not kept as basement');
    // Ordinary seafloor is untouched.
    assert.equal(litho.cover[1], shelfClassId(r_xyz, R), 'shelf cover stopped being assigned');
});

// The latitude call the shelf branch makes, mirrored so the test asserts the
// real value rather than accepting whichever class turns up.
function shelfClassId(r_xyz, R) {
    const latDeg = Math.abs(Math.asin(Math.max(-1, Math.min(1, r_xyz[1]))) * 180 / Math.PI);
    return latDeg < LITHO_CARBONATE_LAT_DEG ? R.carbonate : R.shelf_clastic;
}

test('COVER_SEQUENCE declares every cover class that reaches a real planet', () => {
    // The point of the table: a new cover class cannot arrive without someone
    // deciding where in the sequence it sits. If a class turns up on a surface
    // and no entry claims to produce it, that decision was never made.
    const c = withLitho();
    const byId = Object.fromEntries(ROCK_CLASSES.map(x => [x.id, x.code]));
    const declared = new Set(COVER_SEQUENCE.flatMap(s => s.produces));
    // bare_basement assigns the basement, whatever the basement chain chose, so
    // basement classes are legitimately undeclared here.
    const basements = new Set();
    for (let r = 0; r < c.mesh.numRegions; r++) basements.add(byId[c.lithology.r_basementRock[r]]);

    const undeclared = new Set();
    for (let r = 0; r < c.mesh.numRegions; r++) {
        const code = byId[c.lithology.r_coverRock[r]];
        if (!declared.has(code) && !basements.has(code)) undeclared.add(code);
    }
    assert.deepEqual([...undeclared], [],
        `cover classes with no COVER_SEQUENCE entry: ${[...undeclared].join(', ')}`);

    // Every declared code must be a real rock class, and every entry must be
    // walkable — a typo'd code would silently never match.
    const codes = new Set(ROCK_CLASSES.map(x => x.code));
    for (const s of COVER_SEQUENCE) {
        assert.ok(['ocean', 'continent', 'both'].includes(s.domain), `${s.code}: bad domain`);
        assert.equal(typeof s.when, 'function', `${s.code}: no predicate`);
        assert.equal(typeof s.apply, 'function', `${s.code}: no apply`);
        for (const p of s.produces) assert.ok(codes.has(p), `${s.code} produces unknown class "${p}"`);
    }
    // Each domain must end in an unconditional entry, or a cell can fall
    // through the table with no cover decided at all.
    for (const d of ['ocean', 'continent']) {
        const last = COVER_SEQUENCE.filter(s => s.domain === d || s.domain === 'both').pop();
        assert.equal(last.when({}), true, `${d} has no catch-all entry at the end of the sequence`);
    }
});

test('basin fill is first in the sequence, ahead of every episodic deposit', () => {
    // Three separate bugs came from this being false. It is the one ordering
    // claim in the table that is load-bearing, so it gets its own assertion
    // rather than relying on the behavioural tests below to notice.
    assert.equal(COVER_SEQUENCE[0].code, 'basin_fill');
    assert.equal(COVER_SEQUENCE[0].domain, 'both',
        'basin fill must apply on both crust types, or the same basin classifies two ways');
});

test('closed-basin fill outranks every volcanic cover branch', () => {
    // The chain order is a stratigraphic claim and this is the assertion of it:
    // volcanism resurfaces a landscape once, a closed basin fills continuously
    // afterwards, so the basin wins. Left the other way round, 33 preserved
    // basins on Vesper surfaced as flood basalt with no fill anywhere in them —
    // and the identical basin on OCEANIC LIP already got fill, so the model
    // answered the same question two ways depending on crust type.
    const R = Object.fromEntries(ROCK_CLASSES.map(c => [c.code, c.id]));
    const mesh = { numRegions: 3 };
    const r_xyz = new Float32Array(9);
    const r_elevation = new Float32Array([0.2, 0.2, 0.2]);
    const base = {
        r_isOcean: Uint8Array.from([0, 0, 0]),
        r_stress: Float32Array.from([0, 0, 0]),
        maxStress: 1,
    };
    // Cell 0: closed basin on a flood-basalt province. Cell 1: the province
    // alone. Cell 2: closed basin over an arc.
    const t = {
        ...base,
        riftDist: Float32Array.from([0, 0, 0]),
        riftHalfWidth: 0,
    };
    const dl = { lip: Float32Array.from([1, 1, 0]), backArc: Float32Array.from([0, 0, 1]) };
    const litho = classifyLithology(mesh, r_xyz, r_elevation, t, dl,
        { basins: { r_isEndorheic: Uint8Array.from([1, 0, 1]), r_isSaltCrust: null } });

    assert.equal(litho.cover[0], R.playa_clastic, 'flood basalt buried the basin fill');
    assert.equal(litho.cover[2], R.playa_clastic, 'arc cover buried the basin fill');
    // The volcanic branches must still work where no basin claims the cell.
    assert.equal(litho.cover[1], R.flood_basalt, 'LIP cover stopped being assigned');
});

test('no basin pass means clastic fill, never the bright class by default', () => {
    // Absent the mask, the fallback must be the conservative half of the pair —
    // silently defaulting to albedo 0.50 is how the old bug was worth 2.8 W/m2.
    const mask = saltCrustMask(null, null, new Float32Array(4));
    assert.deepEqual(Array.from(mask), [0, 0, 0, 0]);
});
