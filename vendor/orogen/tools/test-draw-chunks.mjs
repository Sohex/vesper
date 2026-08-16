// Draw-call chunking: the arithmetic behind the viewer's surface split.
//
// The failure this guards against is silent — WebGL logs a warning and draws
// nothing when a single drawArrays exceeds the vertex-id ceiling — so the cap
// is asserted here rather than discovered by looking at a blank canvas.

import test from 'node:test';
import assert from 'node:assert/strict';

import { planChunks, MAX_VERTS_PER_DRAW } from '../js/draw-chunks.js';

/** numSides for a dual mesh of n regions: 3 per triangle, ~2n triangles. */
const sidesFor = (regions) => 6 * regions;

test('no chunk exceeds the per-draw vertex ceiling', () => {
    for (const regions of [2_000, 60_000, 250_000, 1_000_000, 2_500_000, 10_000_000]) {
        const tris = sidesFor(regions);
        for (const c of planChunks(tris)) {
            assert.ok(c.count * 3 <= MAX_VERTS_PER_DRAW,
                `${regions} regions: a chunk asks for ${c.count * 3} vertices`);
        }
    }
});

test('chunks tile the triangle range exactly — no gap, no overlap', () => {
    for (const tris of [0, 1, 999, 2_000_000, 14_981_994]) {
        const chunks = planChunks(tris);
        let cursor = 0;
        for (const c of chunks) {
            assert.equal(c.start, cursor, `gap or overlap at ${c.start}`);
            assert.ok(c.count > 0);
            cursor += c.count;
        }
        assert.equal(cursor, tris, `chunks cover ${cursor} of ${tris} triangles`);
    }
});

test('the reported planet is chunked, and would not have been drawable whole', () => {
    // 2,500,001 regions — the export that failed to render.
    const tris = 14_981_994;
    assert.ok(tris * 3 > 30_000_000, 'fixture no longer reproduces the overflow');
    const chunks = planChunks(tris);
    assert.ok(chunks.length > 1, 'this planet must be split across draws');
    assert.equal(chunks.reduce((n, c) => n + c.count, 0), tris);
});

test('a small planet stays a single draw', () => {
    // Chunking must not fragment the common case into needless draw calls.
    const chunks = planChunks(sidesFor(60_000));
    assert.equal(chunks.length, 1);
    assert.equal(chunks[0].start, 0);
});

test('subarray views over one backing buffer alias by global offset', () => {
    // This is what lets every colour-update path keep writing at s*9 with no
    // per-chunk bookkeeping. If it ever stops being true, colours land in the
    // wrong chunk and the globe renders stale.
    const tris = 10;
    const col = new Float32Array(tris * 9);
    const chunks = planChunks(tris, 12);          // 4 triangles per chunk
    assert.equal(chunks.length, 3);
    const views = chunks.map(c => col.subarray(c.start * 9, (c.start + c.count) * 9));

    col[7 * 9 + 4] = 0.5;                          // triangle 7, a middle chunk
    assert.equal(views[1][(7 - 4) * 9 + 4], 0.5);
    assert.equal(views[0].buffer, col.buffer, 'views must share storage, not copy');
    assert.equal(views.reduce((n, v) => n + v.length, 0), col.length);
});
