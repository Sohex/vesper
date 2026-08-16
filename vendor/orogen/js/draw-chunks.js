// Splitting a surface across several draw calls.
//
// A single drawArrays has a vertex-count ceiling. Firefox enforces 30,000,000
// (webgl.max-vert-ids-per-draw) and drivers impose their own; exceeding it does
// not throw, it logs a warning and draws nothing. The globe is one triangle per
// dual-mesh side, so numSides ≈ 6 × numRegions and a 2.5M-region planet is
// ~15M triangles = ~45M vertices in one call — half again over the cap, and the
// entire surface silently disappeared.
//
// This fork generates at those region counts by design (2K to 2.5M), so the
// viewer splits the surface instead of capping the mesh. Kept DOM-free and
// THREE-free so the arithmetic can be tested headlessly.

/** Vertices per draw call. Well under the 30M ceiling, and a round 2M triangles. */
export const MAX_VERTS_PER_DRAW = 6_000_000;

/**
 * Contiguous triangle ranges covering [0, totalTris), each within the vertex cap.
 *
 * Ranges are returned in order and are exact — every triangle appears in
 * exactly one chunk, so a chunk's float range is [start*9, (start+count)*9) and
 * a write at a global offset lands in exactly one chunk's view.
 *
 * @param {number} totalTris
 * @param {number} maxVerts vertices permitted in one draw call
 * @returns {Array<{start: number, count: number}>}
 */
export function planChunks(totalTris, maxVerts = MAX_VERTS_PER_DRAW) {
    const n = Math.max(0, Math.floor(totalTris));
    if (n === 0) return [];
    const perChunk = Math.max(1, Math.floor(maxVerts / 3));
    const out = [];
    for (let start = 0; start < n; start += perChunk) {
        out.push({ start, count: Math.min(perChunk, n - start) });
    }
    return out;
}
