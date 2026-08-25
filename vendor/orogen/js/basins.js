// Closed-basin (endorheic) detection, cataloguing and protection.
//
// WHAT THIS MODULE DOES AND DELIBERATELY DOES NOT DO
//
// Orogen's job here is to stop destroying depression geometry and to hand
// downstream everything needed to reason about it. It does NOT decide which
// basins are endorheic. Whether a closed basin holds a terminal lake, brims
// over and becomes through-flowing, or sits dry as a salt pan is a water
// balance — catchment precipitation against evaporation over the lake surface
// — and climate belongs to ExoPlaSim, not here. Three separate things, kept
// separate on purpose:
//
//   1. Potential basin geometry  — every closed depression the terrain has.
//      Measured here, exported in full (hypsometry, spill level, catchment).
//   2. Preserved endorheic terrain — the subset whose rims the drainage
//      conditioning is forbidden to breach. Selection is a signal-vs-noise
//      question (is this a resolvable landform?), not a hydrological one.
//   3. Actual lake water level — never decided here. Supplied externally, or
//      left unset. A preserved basin with no water level is a dry closed
//      basin, which is a perfectly good answer.
//
// Everything is deterministic: no Math.random, no Set-iteration-order
// dependence, and every ordering is by region index or an explicit sort.

import { MinHeap } from './min-heap.js';
import { avgCellAreaKm2, avgEdgeKm, uniformCellArea } from './geometry.js';
import { scaledHeightKm } from './color-map.js';
import {
    BASIN_MIN_DEPTH_KM, BASIN_MIN_AREA_KM2, BASIN_MIN_CELLS,
    BASIN_MAX_NEST_DEPTH, BASIN_HYPSOMETRY_LEVELS, BASIN_DIVIDE_RING_KM,
} from './terrain-config.js';

const EPS = 1e-7;

/**
 * Physical height of a point inside a basin, km.
 *
 * Always the land branch: a closed-basin floor below sea level is dry ground,
 * not seabed, and the bathymetric branch reads a -560 m salt pan as a -5.6 km
 * trench. Everything this module labels `...Km` goes through here, so the suffix
 * means what it says.
 *
 * `reliefScale` is the planet's 1/g relief scaling, and it must be applied here
 * for that last sentence to be true. It was not, and the omission was invisible
 * on Earth (scale exactly 1) and only 4% on a 1.04-g planet, but it put the
 * catalogue's depths and volumes on a different vertical scale from the
 * exported `elevation_km` — 31% apart at 1.31 g — while both were published in
 * the same manifest under the same units. A water balance that takes its
 * levels from one and its surface from the other then compares two different
 * kilometres. Positive heights scale and negative ones do not, matching
 * `elevation_km` exactly; the two conversions must stay identical, so if one
 * changes, change both.
 */
const heightKm = (modelElev, reliefScale = 1) =>
    scaledHeightKm(modelElev, reliefScale, true);

/**
 * Sort key for "largest basin first", in km³ at unit relief. See the note in
 * measureDepression: enumeration order, and therefore basin_index, must not
 * move when gravity does. Falls back to volumeKm3 for basin records built by
 * hand (tests) rather than by measureDepression.
 */
const orderingVolume = (b) => b.orderingVolumeKm3 ?? b.volumeKm3;

// ─────────────────────────────────────────────────────────────────────────
//  Ocean topology
// ─────────────────────────────────────────────────────────────────────────

/**
 * Mark the largest connected ocean component as "open ocean". Inland seas are
 * deliberately excluded: a below-sea-level body with no connection to the
 * world ocean is a candidate endorheic basin, not a drainage outlet.
 *
 * Identical to the component scan inside priorityFloodCarve, so both agree on
 * what counts as the sea.
 */
export function findOpenOcean(mesh, r_isOcean) {
    const { numRegions, adjOffset, adjList } = mesh;
    const label = new Int32Array(numRegions).fill(-1);
    const sizes = [];
    for (let r = 0; r < numRegions; r++) {
        if (!r_isOcean[r] || label[r] >= 0) continue;
        const id = sizes.length;
        let size = 0;
        const queue = [r];
        label[r] = id;
        while (queue.length > 0) {
            const cur = queue.pop();
            size++;
            for (let i = adjOffset[cur], iEnd = adjOffset[cur + 1]; i < iEnd; i++) {
                const nb = adjList[i];
                if (r_isOcean[nb] && label[nb] < 0) { label[nb] = id; queue.push(nb); }
            }
        }
        sizes.push(size);
    }
    let main = 0;
    for (let i = 1; i < sizes.length; i++) if (sizes[i] > sizes[main]) main = i;
    const isOpenOcean = new Uint8Array(numRegions);
    if (sizes.length > 0) {
        for (let r = 0; r < numRegions; r++) {
            if (r_isOcean[r] && label[r] === main) isOpenOcean[r] = 1;
        }
    }
    return { isOpenOcean, oceanLabel: label, componentSizes: sizes, mainOceanLabel: main };
}

// ─────────────────────────────────────────────────────────────────────────
//  Fill surface
// ─────────────────────────────────────────────────────────────────────────

/**
 * Barnes et al. priority-flood, computing for each cell the water level at
 * which it becomes submerged. NON-MUTATING — r_elevation is only read.
 *
 * fill[r] === elevation[r] means the cell drains freely to an outlet.
 * fill[r] >  elevation[r] means the cell sits inside a depression, and fill[r]
 * is that depression's spill elevation. No epsilon is added while filling, so
 * every cell of one depression carries a bit-identical fill value — which is
 * what lets components be grouped by exact equality below.
 *
 * @param outlets Uint8Array marking cells that act as drainage outlets
 *                (open ocean, plus any preserved basin sinks).
 */
export function computeFillSurface(mesh, r_elevation, outlets, scope = null, scratch = null, receivers = null) {
    const { numRegions, adjOffset, adjList } = mesh;
    const fill = scratch ? scratch.fill : new Float32Array(numRegions);
    const visited = scratch ? scratch.visited : new Uint8Array(numRegions);
    const heap = new MinHeap(fill);

    // Cells to consider. Passing an explicit list keeps nested detection
    // proportional to the basin rather than to the whole planet.
    const cells = scope ? scope.cells : null;
    const inScope = scope ? scope.mask : null;
    const n = cells ? cells.length : numRegions;
    const at = (i) => (cells ? cells[i] : i);

    for (let i = 0; i < n; i++) {
        const r = at(i);
        visited[r] = 0;
        if (outlets[r]) {
            fill[r] = r_elevation[r];
            visited[r] = 1;
            if (receivers) receivers[r] = -1;      // an outlet is where water stops
            heap.push(r);
        }
    }
    while (heap.size > 0) {
        const r = heap.pop();
        const level = fill[r];
        for (let i = adjOffset[r], iEnd = adjOffset[r + 1]; i < iEnd; i++) {
            const nb = adjList[i];
            if (inScope && !inScope[nb]) continue;   // outside the scope is a wall
            if (visited[nb]) continue;
            visited[nb] = 1;
            fill[nb] = r_elevation[nb] > level ? r_elevation[nb] : level;
            // The cell that first reached nb is the one nb drains to. Collected
            // over the whole flood this is a spanning tree rooted at the
            // outlets, so it routes across filled flats where steepest descent
            // has no gradient to follow at all.
            if (receivers) receivers[nb] = r;
            heap.push(nb);
        }
    }
    // Cells unreachable from any outlet (no outlet at all, or an isolated
    // island of mesh) keep their own elevation — they are not depressions.
    for (let i = 0; i < n; i++) {
        const r = at(i);
        if (!visited[r]) fill[r] = r_elevation[r];
    }
    return fill;
}

// ─────────────────────────────────────────────────────────────────────────
//  Stable identity
// ─────────────────────────────────────────────────────────────────────────

/** FNV-1a 32-bit over a sorted Int32 sequence. */
function fnv1a32(sortedMembers) {
    let h = 0x811c9dc5;
    for (let i = 0; i < sortedMembers.length; i++) {
        let v = sortedMembers[i];
        for (let b = 0; b < 4; b++) {
            h ^= (v >>> (b * 8)) & 0xff;
            h = Math.imul(h, 0x01000193) >>> 0;
        }
    }
    return h >>> 0;
}

/**
 * Deterministic basin identity: sink region plus a hash of the sorted member
 * set. Stable across runs of the same seed and parameters.
 *
 * NOT stable across region counts — region indices are a property of the mesh,
 * so the same planet at 100k and 250k regions gives a basin different IDs.
 * Match basins across resolutions by centroid and area, not by ID.
 */
export function basinId(sink, sortedMembers) {
    return `orogen-basin-r${sink.toString(36)}-${fnv1a32(sortedMembers).toString(36)}`;
}

// ─────────────────────────────────────────────────────────────────────────
//  Depression detection
// ─────────────────────────────────────────────────────────────────────────

/**
 * Group submerged cells into depressions and measure each one.
 *
 * Cells belong to the same depression when they are adjacent AND share a fill
 * value exactly — adjacent depressions that spill at different levels stay
 * distinct. Nested sub-depressions are merged into their parent at this level;
 * detectBasins() recovers the nesting by recursing.
 *
 * @param restrict optional Uint8Array; when given, only cells marked 1 are
 *                 considered (used when recursing inside a parent basin).
 */
function groupDepressions(mesh, r_elevation, fill, scope, component) {
    const { adjOffset, adjList } = mesh;
    const groups = [];

    const cells = scope ? scope.cells : null;
    const inScope = scope ? scope.mask : null;
    const n = cells ? cells.length : mesh.numRegions;
    const at = (i) => (cells ? cells[i] : i);

    for (let i = 0; i < n; i++) component[at(i)] = -1;

    const isDepression = (r) => fill[r] > r_elevation[r] + EPS;

    for (let i = 0; i < n; i++) {
        const seed = at(i);
        if (component[seed] >= 0 || !isDepression(seed)) continue;
        const level = fill[seed];
        const id = groups.length;
        const members = [];
        const queue = [seed];
        component[seed] = id;
        while (queue.length > 0) {
            const r = queue.pop();
            members.push(r);
            for (let j = adjOffset[r], jEnd = adjOffset[r + 1]; j < jEnd; j++) {
                const nb = adjList[j];
                if (inScope && !inScope[nb]) continue;
                if (component[nb] >= 0 || !isDepression(nb)) continue;
                if (fill[nb] !== level) continue;
                component[nb] = id;
                queue.push(nb);
            }
        }
        members.sort((a, b) => a - b);
        groups.push({ level, members });
    }
    return groups;
}

/** Measure one depression: sink, depth, area, volume, spill point, outlet. */
function measureDepression(mesh, r_elevation, group, component, cellArea, componentId, reliefScale = 1) {
    const { adjOffset, adjList } = mesh;
    const { level, members } = group;

    let sink = members[0];
    let sinkElev = r_elevation[sink];
    let area = 0;
    let volume = 0;
    for (let i = 0; i < members.length; i++) {
        const r = members[i];
        const e = r_elevation[r];
        if (e < sinkElev) { sinkElev = e; sink = r; }
        const a = cellArea[r];
        area += a;
        const d = level - e;
        if (d > 0) volume += d * a;
    }

    // Spill point: the member/non-member pair whose higher side is lowest —
    // the saddle the depression overtops first. Ties break on region index so
    // the choice is deterministic.
    let spillFrom = -1, spillTo = -1, spillSaddle = Infinity;
    for (let i = 0; i < members.length; i++) {
        const r = members[i];
        for (let j = adjOffset[r], jEnd = adjOffset[r + 1]; j < jEnd; j++) {
            const nb = adjList[j];
            if (component[nb] === componentId) continue;
            const saddle = Math.max(r_elevation[r], r_elevation[nb]);
            if (saddle < spillSaddle ||
                (saddle === spillSaddle && (r < spillFrom || (r === spillFrom && nb < spillTo)))) {
                spillSaddle = saddle;
                spillFrom = r;
                spillTo = nb;
            }
        }
    }

    // Physical volume: integrate the height difference, not the model-parameter
    // difference. The conversion is nonlinear, so scaling a model-unit volume
    // afterwards would be wrong.
    const spillH = heightKm(level, reliefScale), sinkH = heightKm(sinkElev, reliefScale);
    let volumeKm3 = 0;
    for (let i = 0; i < members.length; i++) {
        const d = spillH - heightKm(r_elevation[members[i]], reliefScale);
        if (d > 0) volumeKm3 += d * cellArea[members[i]];
    }

    // Enumeration order must not depend on gravity. selectBasins sorts "largest
    // basin first", and that order assigns basin_index and the order of the
    // preserved list. Sorting on the g-scaled volume let a gravity change
    // reorder two basins of near-equal size — the scaling is not a single
    // factor, because heights below sea level do not scale — which moved
    // basin_index and the basinsPreserved hash on a planet whose basins are
    // identical. Ids are content-derived and stayed put, so the reshuffle was
    // invisible to anything joining on them and visible to anything joining on
    // the index. Same volume measured at unit relief: identical at Earth
    // gravity, and identical to the ordering every build before the 1/g fix
    // used, so that fix changes published numbers and nothing structural.
    let orderingVolumeKm3 = volumeKm3;
    if (reliefScale !== 1) {
        const spillH1 = heightKm(level);
        orderingVolumeKm3 = 0;
        for (let i = 0; i < members.length; i++) {
            const d = spillH1 - heightKm(r_elevation[members[i]]);
            if (d > 0) orderingVolumeKm3 += d * cellArea[members[i]];
        }
    }

    return {
        sink,
        // Model-parameter values, unsuffixed so the name does not lie.
        sinkElevation: sinkElev,
        spillElevation: level,
        depth: level - sinkElev,
        // Physical values. These are the ones a water balance wants.
        sinkElevationKm: sinkH,
        spillElevationKm: spillH,
        depthKm: spillH - sinkH,
        areaKm2: area,
        volumeKm3,
        // Internal, deliberately not in the manifest whitelist: a sort key, not
        // a measurement. Equals volumeKm3 at Earth gravity.
        orderingVolumeKm3,
        volumeModelUnits: volume,
        cellCount: members.length,
        spillFrom,
        spillTo,
        spillSaddle: spillSaddle === Infinity ? null : spillSaddle,
        spillSaddleKm: spillSaddle === Infinity ? null : heightKm(spillSaddle, reliefScale),
        members,
    };
}

/**
 * Hypsometry: how area and volume grow with water level, from the sink up to
 * the spill point. This is what a downstream water balance integrates against
 * — it is the basin's actual shape, independent of whether any water is in it.
 */
export function basinHypsometry(r_elevation, members, cellArea, sinkElev, spillElev, levels, reliefScale = 1) {
    const n = Math.max(2, levels | 0);
    const out = [];
    const span = spillElev - sinkElev;
    const sinkH = heightKm(sinkElev, reliefScale);
    for (let i = 0; i < n; i++) {
        // Sample evenly in the MODEL parameter (that is where the terrain lives),
        // but report and integrate in physical km — the conversion is nonlinear,
        // so a volume computed in model units cannot be rescaled afterwards.
        const level = sinkElev + span * (i / (n - 1));
        const levelH = heightKm(level, reliefScale);
        let area = 0, volumeKm3 = 0;
        for (let k = 0; k < members.length; k++) {
            const e = r_elevation[members[k]];
            if (e <= level) {
                const a = cellArea[members[k]];
                area += a;
                volumeKm3 += (levelH - heightKm(e, reliefScale)) * a;
            }
        }
        out.push({
            level: level,
            levelKm: levelH,
            depthKm: levelH - sinkH,
            areaKm2: area,
            volumeKm3,
        });
    }
    return out;
}

/**
 * Full depression catalogue for a terrain, with nesting.
 *
 * Top-level basins are the depressions of the terrain as a whole. Each is then
 * re-flooded from its own sink to find the sub-depressions nested inside it,
 * down to BASIN_MAX_NEST_DEPTH. The result is a forest: every basin knows its
 * parent and the basin (or the ocean) it spills into.
 *
 * NON-MUTATING.
 */
export function detectBasins(mesh, r_elevation, r_isOcean, opts = {}) {
    const { numRegions } = mesh;
    const radiusKm = opts.radiusKm;
    const cellArea = opts.cellArea || uniformCellArea(numRegions, radiusKm);
    const maxNest = opts.maxNestDepth ?? BASIN_MAX_NEST_DEPTH;
    // Threaded, never module state: several pipelines run in one process.
    const reliefScale = opts.reliefScale ?? 1;

    const { isOpenOcean } = opts.isOpenOcean
        ? { isOpenOcean: opts.isOpenOcean }
        : findOpenOcean(mesh, r_isOcean);

    // Scratch buffers allocated once and reused at every level of the
    // recursion. Allocating per basin was O(basins x numRegions) in both time
    // and garbage, which on a real 500k-region planet with thousands of
    // depressions dominated the entire pipeline.
    const scratch = { fill: new Float32Array(numRegions), visited: new Uint8Array(numRegions) };
    const component = new Int32Array(numRegions);
    const scopeMask = new Uint8Array(numRegions);
    const outlets = new Uint8Array(numRegions);

    const basins = [];

    /**
     * @param scope null for the whole planet, else {cells, mask} for one basin.
     * Work at each level is proportional to the cells in scope, so the total
     * across a level is O(numRegions log numRegions) however many basins there are.
     */
    const recurse = (scope, parentIndex, depth) => {
        const fill = computeFillSurface(mesh, r_elevation, outlets, scope, scratch);
        const groups = groupDepressions(mesh, r_elevation, fill, scope, component);
        const produced = [];
        for (let g = 0; g < groups.length; g++) {
            const m = measureDepression(mesh, r_elevation, groups[g], component, cellArea, g, reliefScale);
            const index = basins.length;
            const record = {
                index,
                id: basinId(m.sink, m.members),
                parentIndex,
                nestDepth: depth,
                childIndices: [],
                ...m,
                // Hypsometry is deliberately NOT computed here. Only preserved
                // basins export it, and a large planet detects tens of thousands
                // of depressions — computing 24 level samples for every one of
                // them and discarding almost all is pure waste. attachHypsometry()
                // fills it in after selection.
            };
            basins.push(record);
            if (parentIndex >= 0) basins[parentIndex].childIndices.push(index);
            produced.push(record);
        }
        if (depth >= maxNest) return;
        for (const record of produced) {
            // Re-flood inside this basin only, with its own sink as the sole
            // outlet: a sub-depression is terrain that cannot drain to the sink.
            // Everything outside the member set is a wall, not an outlet.
            for (const r of record.members) scopeMask[r] = 1;
            outlets[record.sink] = 1;
            recurse({ cells: record.members, mask: scopeMask }, record.index, depth + 1);
            outlets[record.sink] = 0;
            for (const r of record.members) scopeMask[r] = 0;
        }
    };

    for (let r = 0; r < numRegions; r++) outlets[r] = isOpenOcean[r];
    recurse(null, -1, 0);
    outlets.fill(0);

    // Resolve what each basin spills into, now that every basin is catalogued.
    const owner = new Int32Array(numRegions).fill(-1);
    for (const b of basins) {
        // Deeper (more nested) basins claim their cells last, so `owner` ends
        // up holding the innermost basin containing each cell.
        for (const r of b.members) {
            if (owner[r] < 0 || basins[owner[r]].nestDepth < b.nestDepth) owner[r] = b.index;
        }
    }
    for (const b of basins) {
        if (b.spillTo < 0) { b.spillsInto = null; continue; }
        if (isOpenOcean[b.spillTo]) { b.spillsInto = 'ocean'; continue; }
        const target = owner[b.spillTo];
        b.spillsInto = target >= 0 && target !== b.index ? basins[target].id
                     : (b.parentIndex >= 0 ? basins[b.parentIndex].id : 'ocean');
    }

    return basins;
}

/** Fill in hypsometry for a chosen set of basins, after selection. */
export function attachHypsometry(selected, r_elevation, cellArea, levels = BASIN_HYPSOMETRY_LEVELS, reliefScale = 1) {
    for (const b of selected) {
        b.hypsometry = basinHypsometry(r_elevation, b.members, cellArea,
            b.sinkElevation, b.spillElevation, levels, reliefScale);
    }
    return selected;
}

// ─────────────────────────────────────────────────────────────────────────
//  Selection
// ─────────────────────────────────────────────────────────────────────────

/**
 * Choose which depressions the drainage conditioning must not breach.
 *
 * This is explicitly NOT a judgement about endorheism — see the module header.
 * It answers one question: is this depression a landform the mesh actually
 * resolves, or is it noise? Three floors, and a depression must clear all of
 * them. They are NOT in the same currency, and that decides which one governs a
 * given build:
 *
 *   area       — BASIN_MIN_AREA_KM2 across, in real km². A fixed physical size,
 *                so the same planet preserves the same basins at any region
 *                count once the mesh can hold them. This is the floor that
 *                makes the feature scale-invariant.
 *   cells      — and at least BASIN_MIN_CELLS cells, because a depression
 *                spanning three cells is a mesh artifact whatever its area
 *                works out to. Binds only while the mesh is coarse enough that
 *                BASIN_MIN_CELLS cells cover more ground than
 *                BASIN_MIN_AREA_KM2.
 *   depth      — and BASIN_MIN_DEPTH_KM below its spill, compared against the
 *                depression's PHYSICAL depth, `b.depthKm`. The constant, the
 *                `resolution.minDepthKm` manifest key and the published
 *                `selectionCriteria.minDepthKm` all name a physical depth, and
 *                this is where one is enforced, so the declared floor means the
 *                same drop wherever the depression sits. Comparing against
 *                `b.depth`, the depression's extent in the model's
 *                DIMENSIONLESS elevation parameter, would not: that curve is
 *                quartic above sea level, linear below it, and saturates at
 *                model elevation 1, so a single model-unit threshold admits
 *                depressions a few metres deep near sea level, demands hundreds
 *                of metres on high ground, and admits depressions lying wholly
 *                above the saturation whose physical depth — and therefore
 *                whose lake capacity — is exactly zero. Rejecting mesh-scale
 *                terrain noise is the area and cell floors' job: a depression
 *                clearing BASIN_MIN_AREA_KM2 and BASIN_MIN_CELLS is terrain the
 *                mesh resolves whatever currency its depth is read in.
 *                `basinResolutionContext` publishes `minDepthComparedIn` so a
 *                consumer of the manifest never has to infer the currency, and
 *                `hydrography/scripts/catalogue_floor.py` measures what the
 *                floor costs on a given build.
 *
 * Explicit IDs bypass every floor: naming a basin means you want it.
 */
/**
 * Parse a drainage-hypothesis file.
 *
 * One basin id per line, optionally followed by a retain fraction in [0,1]:
 *
 *     orogen-basin-r1a2-xyz          # fully preserved
 *     orogen-basin-r9zz-abc  0.35    # rim may lose 65% of the basin's relief
 *
 * `#` starts a comment. A JSON array of ids, or of {id, retain} objects, is
 * also accepted — the water-balance step that produces these lists is likelier
 * to emit JSON than a text file.
 */
export function parseBasinList(text) {
    const trimmed = String(text).trim();
    if (trimmed.startsWith('[') || trimmed.startsWith('{')) {
        const parsed = JSON.parse(trimmed);
        const arr = Array.isArray(parsed) ? parsed : (parsed.basins || parsed.ids || []);
        return arr.map(x => (typeof x === 'string'
            ? { id: x, retain: 1 }
            : { id: x.id, retain: x.retain === undefined ? 1 : Number(x.retain) }))
            .filter(x => x.id);
    }
    const out = [];
    for (const rawLine of trimmed.split(/\r?\n/)) {
        const line = rawLine.replace(/#.*$/, '').trim();
        if (!line) continue;
        const parts = line.split(/[\s,]+/);
        const id = parts[0];
        const retain = parts.length > 1 ? Number(parts[1]) : 1;
        if (!Number.isFinite(retain) || retain < 0 || retain > 1) {
            throw new Error(`Basin list: retain fraction for "${id}" must be in [0,1], got "${parts[1]}"`);
        }
        out.push({ id, retain });
    }
    return out;
}

export function selectBasins(basins, opts = {}) {
    const minDepthKm = opts.minDepthKm ?? BASIN_MIN_DEPTH_KM;
    const minAreaKm2 = opts.minAreaKm2 ?? BASIN_MIN_AREA_KM2;
    const minCells   = opts.minCells   ?? BASIN_MIN_CELLS;
    const explicitIds = new Set(opts.explicitIds || []);
    const enabled = opts.enabled !== false;

    // A drainage hypothesis from downstream. `preserveList` replaces the
    // threshold selection outright — the preserved set becomes exactly what was
    // asked for, which is what makes an iterative water-balance loop converge.
    // `carveList` is subtractive: threshold selection minus these.
    const preserveList = opts.preserveList || null;
    const carveList = opts.carveList || null;
    const retainById = new Map();
    for (const e of (preserveList || [])) retainById.set(e.id, e.retain);
    const carveIds = new Map((carveList || []).map(e => [e.id, e.retain]));

    const includeNested = opts.includeNested === true;
    const warnings = [];
    const selected = [];
    const byId = new Map(basins.map(b => [b.id, b]));

    for (const id of explicitIds) {
        if (!byId.has(id)) warnings.push(`No basin with id "${id}" in this terrain — ignored.`);
    }

    // Threshold-qualifying basins, outermost first, so the nesting filter below
    // can ask "is an ancestor already in?".
    const qualifies = (b) => enabled
        // `b.depthKm`, not `b.depth`: the floor is declared, named and published
        // as a physical depth, so it is compared as one. The docstring above has
        // the argument and says what the model-unit comparison cost.
        && b.depthKm >= minDepthKm
        && b.areaKm2 >= minAreaKm2
        && b.cellCount >= minCells;

    // An explicit preserve list short-circuits the floors entirely.
    if (preserveList) {
        const byId = new Map(basins.map(b => [b.id, b]));
        let carvedByList = 0;
        for (const e of preserveList) {
            const b = byId.get(e.id);
            if (!b) { warnings.push(`Preserve list names "${e.id}", which is not in this catalogue — ignored.`); continue; }
            // retain 0 is a CARVE verdict, not a preservation with a full-depth
            // notch. Selecting it anyway would make the basin a terminal drainage
            // root — water routes inward to a sink whose rim is unprotected — which
            // is a state the terrain cannot mean. A complete verdict listing every
            // basin (carved ones at 0 rather than omitted) is the natural thing for
            // downstream to send, so it must mean here what it says there. The
            // carve-list path has always dropped these; this one did not.
            if (!(e.retain > 0)) { carvedByList++; continue; }
            selected.push({ ...b, selectedBy: 'preserve-list', retain: e.retain });
        }
        selected.sort((a, b) => orderingVolume(b) - orderingVolume(a) || a.sink - b.sink);
        return { selected, warnings, source: 'preserve-list', carvedByList };
    }

    const chosen = new Set();
    const ordered = basins.slice().sort((a, b) => a.nestDepth - b.nestDepth || a.index - b.index);
    for (const b of ordered) {
        const explicit = explicitIds.has(b.id);
        if (!explicit && !qualifies(b)) continue;
        // Carve list: drop it, or keep it with a reduced retain fraction.
        if (carveIds.has(b.id)) {
            const retain = carveIds.get(b.id);
            if (!(retain > 0)) continue;
            chosen.add(b.index);
            selected.push({ ...b, selectedBy: 'threshold-carved', retain });
            continue;
        }
        if (!explicit && !includeNested) {
            // A sub-basin inside an already-preserved basin adds a divide
            // *within* a preserved basin, which is almost never what a size
            // threshold was meant to express. Explicit ids still get through.
            let p = b.parentIndex, skip = false;
            while (p >= 0) {
                if (chosen.has(p)) { skip = true; break; }
                p = basins[p].parentIndex;
            }
            if (skip) continue;
        }
        chosen.add(b.index);
        selected.push({ ...b, selectedBy: explicit ? 'explicit' : 'threshold' });
    }

    // Nested selection can still happen via explicit ids or includeNested. It
    // is not an error, but the inner rim becomes a divide inside a preserved
    // basin, so say so rather than doing it silently.
    const selectedIndices = new Set(selected.map(b => b.index));
    for (const b of selected) {
        let p = b.parentIndex;
        while (p >= 0) {
            if (selectedIndices.has(p)) {
                warnings.push(
                    `Basin ${b.id} (nest depth ${b.nestDepth}) is nested inside also-selected ` +
                    `${basins[p].id}. Both rims will be protected; the inner one is a divide ` +
                    `inside a preserved basin.`);
                break;
            }
            p = basins[p].parentIndex;
        }
    }

    for (const id of carveIds.keys()) {
        if (!basins.some(b => b.id === id)) {
            warnings.push(`Carve list names "${id}", which is not in this catalogue — ignored.`);
        }
    }

    selected.sort((a, b) => orderingVolume(b) - orderingVolume(a) || a.sink - b.sink);
    return { selected, warnings, source: carveList ? 'threshold-minus-carve-list' : 'threshold' };
}

// ─────────────────────────────────────────────────────────────────────────
//  Protection
// ─────────────────────────────────────────────────────────────────────────

/**
 * Build the masks the drainage conditioning consults.
 *
 *   terminalRoot — each preserved basin's sink. priorityFloodCarve seeds its
 *                  flood from these as well as from the ocean, so basin cells
 *                  drain inward to the sink instead of being filled and carved
 *                  outward toward the sea.
 *   noLower      — the basin divide: member cells on the rim plus a ring of
 *                  their outside neighbours. The carve passes exist precisely
 *                  to cut outlets through depressions, and left alone they take
 *                  ~900 m out of a rim. They must not touch these.
 *   carveAllowance — how far below its running elevation a divide cell may be
 *                  carved: 0 for a fully preserved rim, (1-retain)·depth for a
 *                  partially preserved one.
 *   basinIndex   — which preserved basin owns each cell, or -1.
 *
 * Reads no elevations. The floors are relative, so this can be built before
 * erosion and still mean the same thing during it.
 */
export function buildBasinProtection(mesh, selected, opts = {}) {
    const { numRegions, adjOffset, adjList } = mesh;
    const r_elevation = opts.r_elevation || null;
    // Band width from a physical distance, not a hop count — see
    // BASIN_DIVIDE_RING_KM. At least one hop, so the band never vanishes.
    const ringKm = opts.ringKm ?? BASIN_DIVIDE_RING_KM;
    const ringCells = opts.ringCells
        ?? Math.max(1, Math.round(ringKm / avgEdgeKm(numRegions, opts.radiusKm)));
    const terminalRoot = new Uint8Array(numRegions);
    const noLower = new Uint8Array(numRegions);
    // How far below its CURRENT elevation each divide cell may be carved. Zero
    // for a fully preserved rim; a partially preserved one gets a slice of its
    // basin's relief, so the carve can incise part-way and leave a
    // through-flowing valley with a residual lake instead of a cut trench.
    //
    // An allowance, not an absolute floor. Protection is built before erosion
    // and the carve runs during it, so an absolute floor pinned here is a
    // pre-erosion height: where erosion had RAISED a divide, the carve was then
    // free to take it back down, and the invariant guard — which uses the same
    // array as its baseline — could not see it. It also moved terrain that
    // nothing had asked to move. Relative to the running surface, retain = 1 is
    // bit-identical to no allowance at all, which is what it should be.
    const carveAllowance = new Float32Array(numRegions);
    const basinIndex = new Int32Array(numRegions).fill(-1);
    const isMember = new Uint8Array(numRegions);
    // retain[i] in [0,1] per selected basin: 1 preserves the rim outright, 0 is
    // equivalent to not selecting it at all.
    const retainOf = (i) => {
        const v = selected[i].retain;
        return v === undefined || v === null ? 1 : Math.max(0, Math.min(1, v));
    };

    for (let i = 0; i < selected.length; i++) {
        for (const r of selected[i].members) {
            isMember[r] = 1;
            // Innermost selection wins when preserved basins nest.
            if (basinIndex[r] < 0 || selected[i].nestDepth > selected[basinIndex[r]].nestDepth) {
                basinIndex[r] = i;
            }
        }
    }
    for (let i = 0; i < selected.length; i++) terminalRoot[selected[i].sink] = 1;

    // Seed the divide with member cells that touch a non-member.
    const frontier = [];
    for (let r = 0; r < numRegions; r++) {
        if (!isMember[r]) continue;
        for (let i = adjOffset[r], iEnd = adjOffset[r + 1]; i < iEnd; i++) {
            if (!isMember[adjList[i]]) { noLower[r] = 1; frontier.push(r); break; }
        }
    }
    // Grow it outward by ringCells hops. The carve kernel has a radius of
    // several cells, so a one-cell divide would just be cut around.
    let current = frontier;
    for (let hop = 0; hop < ringCells; hop++) {
        const next = [];
        for (const r of current) {
            for (let i = adjOffset[r], iEnd = adjOffset[r + 1]; i < iEnd; i++) {
                const nb = adjList[i];
                if (!noLower[nb]) { noLower[nb] = 1; next.push(nb); }
            }
        }
        current = next;
    }
    // A sink is never a divide — protecting it would block legitimate erosion
    // of the basin floor, which we explicitly leave enabled.
    for (let i = 0; i < selected.length; i++) noLower[selected[i].sink] = 0;

    // Allowance is a fraction of the owning basin's natural depth, so
    // "retain 0.8" means the rim may lose a fifth of the basin's relief — a
    // proportional incision, not a breach.
    for (let r = 0; r < numRegions; r++) {
        if (!noLower[r]) continue;
        const bi = basinIndex[r];
        const retain = bi >= 0 ? retainOf(bi) : 1;
        if (retain >= 1) continue;                       // allowance stays 0
        if (retain <= 0) { noLower[r] = 0; continue; }   // not protected at all
        const depth = bi >= 0 ? (selected[bi].depth ?? 0) : 0;
        carveAllowance[r] = (1 - retain) * depth;
    }

    return { terminalRoot, noLower, carveAllowance, basinIndex, isMember };
}

/**
 * Cut partial outlets through the rims of basins whose retain fraction is < 1.
 *
 * MUTATES r_elevation.
 *
 * Lowering the protection floor alone does nothing: with a terminal root the
 * flood has no deficit at that basin, so the carve never fires and the floor is
 * never approached. Permission is not incision — the notch has to be cut
 * deliberately.
 *
 * The cut is a notch at the basin's saddle, deepest at the spill point and
 * ramping back up over the divide band, lowering the rim by (1 - retain) of the
 * basin's relief. A basin that only just overflows therefore ends up a
 * through-flowing valley holding a residual lake, rather than either an intact
 * closed basin or a fully trenched one.
 */
export function inciseOutlets(mesh, r_elevation, selected, opts = {}) {
    const { adjOffset, adjList } = mesh;
    const hops = Math.max(1, opts.ringCells
        ?? Math.round(BASIN_DIVIDE_RING_KM / avgEdgeKm(mesh.numRegions, opts.radiusKm)));
    let cut = 0, basinsCut = 0;

    for (const b of selected) {
        const retain = b.retain === undefined || b.retain === null ? 1 : b.retain;
        if (retain >= 1 || b.spillFrom < 0) continue;
        const allowance = (1 - retain) * (b.depth ?? 0);
        if (!(allowance > 0)) continue;
        basinsCut++;

        // BFS out from the saddle; the notch is deepest at the spill point and
        // tapers linearly to nothing at the edge of the divide band.
        const seen = new Map([[b.spillFrom, 0]]);
        const queue = [b.spillFrom];
        if (b.spillTo >= 0) { seen.set(b.spillTo, 0); queue.push(b.spillTo); }
        for (let qi = 0; qi < queue.length; qi++) {
            const r = queue[qi];
            const d = seen.get(r);
            const target = (b.spillElevation ?? r_elevation[r]) - allowance * (1 - d / (hops + 1));
            if (r_elevation[r] > target) { r_elevation[r] = target; cut++; }
            if (d >= hops) continue;
            for (let i = adjOffset[r], iEnd = adjOffset[r + 1]; i < iEnd; i++) {
                const nb = adjList[i];
                if (!seen.has(nb)) { seen.set(nb, d + 1); queue.push(nb); }
            }
        }
    }
    return { basinsCut, cellsLowered: cut };
}

// ─────────────────────────────────────────────────────────────────────────
//  Final-terrain drainage
// ─────────────────────────────────────────────────────────────────────────

/**
 * Steepest-descent drainage on a finished terrain, with preserved basin sinks
 * acting as terminals alongside the ocean.
 *
 * The catchment of the terrain that actually shipped is not the same as the
 * catchment of the pre-conditioning terrain — erosion moves divides. Hydrology
 * downstream wants this one.
 */
export function computeFinalTerrainDrainage(mesh, r_elevation, r_isOcean, terminalRoot, opts = {}) {
    const { numRegions, adjOffset, adjList } = mesh;
    // Only the connected world ocean is a drainage terminal. Using
    // elevation <= 0 here would make a deep endorheic basin's own floor a
    // terminal and shatter its catchment — the exact distinction this whole
    // module exists to preserve.
    const { isOpenOcean } = findOpenOcean(mesh, r_isOcean);

    // Route on the FILLED surface unless told otherwise.
    //
    // Raw steepest descent on a finished terrain strands every cell that drains
    // into a local pit: erosion creates pits after the last conditioning pass,
    // and their catchments resolve to no terminal at all. On a real planet that
    // was 63% of the land, so a catchment measured this way counted less than
    // half its true area — while being documented as the area to integrate
    // precipitation over. Filling first routes through pits to their spill
    // points, which is what a drop of water actually does. Pass {raw:true} for
    // the unrouted measure.
    const drainTo = new Int32Array(numRegions).fill(-1);
    const terminal = new Int32Array(numRegions).fill(-1);
    const isTerminal = (r) => isOpenOcean[r] || (terminalRoot && terminalRoot[r]);

    let routeSurface = r_elevation;
    if (!opts.raw) {
        // Filling alone is not enough: a filled depression is a plateau, and
        // steepest descent has no gradient to follow across it. The flood's own
        // visitation tree does — every cell records the neighbour that first
        // reached it, which is a path to an outlet by construction.
        const outlets = new Uint8Array(numRegions);
        for (let r = 0; r < numRegions; r++) if (isTerminal(r)) outlets[r] = 1;
        const receivers = new Int32Array(numRegions).fill(-1);
        routeSurface = computeFillSurface(mesh, r_elevation, outlets, null, null, receivers);
        for (let r = 0; r < numRegions; r++) {
            drainTo[r] = isTerminal(r) ? -1 : receivers[r];
        }
    }

    const order = [];
    for (let r = 0; r < numRegions; r++) order.push(r);
    order.sort((a, b) => routeSurface[a] - routeSurface[b] || a - b);

    for (let r = 0; r < numRegions; r++) {
        if (isTerminal(r)) { terminal[r] = r; continue; }
        if (!opts.raw) continue;                 // drainTo already set from the flood tree
        let best = -1, bestElev = r_elevation[r];
        for (let i = adjOffset[r], iEnd = adjOffset[r + 1]; i < iEnd; i++) {
            const nb = adjList[i];
            if (r_elevation[nb] < bestElev) { bestElev = r_elevation[nb]; best = nb; }
        }
        drainTo[r] = best;
    }

    // Resolve each cell to its terminal by walking downhill, ascending
    // elevation order so targets are already resolved.
    for (let k = 0; k < order.length; k++) {
        const r = order[k];
        if (terminal[r] >= 0) continue;
        const chain = [];
        let cur = r;
        while (cur >= 0 && terminal[cur] < 0) {
            chain.push(cur);
            if (chain.length > numRegions) break;   // cycle guard
            cur = drainTo[cur];
        }
        const end = cur >= 0 && terminal[cur] >= 0 ? terminal[cur] : -2;  // -2 = endorheic-but-unclaimed / cycle
        for (const c of chain) terminal[c] = end;
    }

    return { drainTo, terminal };
}

/**
 * Catchment of each preserved basin on a given surface: every cell whose
 * routed path ends at that basin's sink. This is the area a water balance
 * integrates precipitation over — usually much larger than the basin floor.
 *
 * ROOT AT THE SINK OF *THIS* SURFACE. `selected[i].sink` is the low point of
 * the pre-conditioning catalogue, and erosion moves the low point inside the
 * member set: on a 2.5M-region planet 878 of 3629 preserved basins ended with a
 * new deepest cell, in every case strictly below the catalogue sink. Rooting
 * the tree at the old cell put `drainage_terminal` on a cell that water flows
 * *out of* on the shipped terrain, and left two published numbers contradicting
 * each other — summing the field by `finalPreserved.sink` gave 58% of land
 * endorheic where `finalCatchment.areaKm2` said 76%, because those 878 roots
 * accumulated nothing at all.
 *
 * So pass `opts.sinks` (one region index per selected basin) whenever the
 * surface is not the one the catalogue was built on. The default is the
 * catalogue sink, which is right for measuring the pre-conditioning surface and
 * wrong for anything downstream of erosion.
 */
export function computeFinalPreservedCatchments(mesh, r_elevation, r_isOcean, selected, cellArea, opts = {}) {
    const { numRegions } = mesh;
    const area = cellArea || uniformCellArea(numRegions);

    const sinks = opts.sinks
        ? Array.from(opts.sinks, (s, i) => (Number.isInteger(s) && s >= 0 ? s : selected[i].sink))
        : selected.map(b => b.sink);
    if (sinks.length !== selected.length) {
        throw new Error(`opts.sinks has ${sinks.length} entries for ${selected.length} basins`);
    }

    // One cell cannot be the terminal of two basins. Nested preserved basins can
    // share a low point; the innermost owns it, matching how buildBasinProtection
    // resolves overlapping membership.
    const sinkToBasin = new Map();
    const sharedRoots = [];
    for (let i = 0; i < sinks.length; i++) {
        const prev = sinkToBasin.get(sinks[i]);
        if (prev === undefined) { sinkToBasin.set(sinks[i], i); continue; }
        sharedRoots.push({ region: sinks[i], ids: [selected[prev].id, selected[i].id] });
        if ((selected[i].nestDepth ?? 0) > (selected[prev].nestDepth ?? 0)) sinkToBasin.set(sinks[i], i);
    }

    // A root outside its own basin is meaningless — it would collect a
    // catchment that has nothing to do with the depression the record names.
    // (Being its own terminal is NOT worth asserting: the roots are exactly the
    // terminal set, so terminal[root] === root by construction. What actually
    // went wrong was rooting at a cell that was never in this set at all, which
    // no invariant inside this function can see. summariseBasins checks that.)
    const terminalRoot = new Uint8Array(numRegions);
    for (let i = 0; i < sinks.length; i++) {
        if (opts.sinks && !selected[i].members.includes(sinks[i])) {
            throw new Error(`Basin ${selected[i].id}: root ${sinks[i]} is not one of its members`);
        }
        terminalRoot[sinks[i]] = 1;
    }

    const { drainTo, terminal } = computeFinalTerrainDrainage(mesh, r_elevation, r_isOcean, terminalRoot);

    const catchmentCells = selected.map(() => []);
    const catchmentAreaKm2 = new Float64Array(selected.length);
    for (let r = 0; r < numRegions; r++) {
        const t = terminal[r];
        if (t < 0) continue;
        const bi = sinkToBasin.get(t);
        if (bi === undefined) continue;
        catchmentCells[bi].push(r);
        catchmentAreaKm2[bi] += area[r];
    }

    const out = selected.map((b, i) => ({
        index: i,
        id: b.id,
        sink: sinks[i],
        catalogueSink: b.sink,
        catchmentCells: catchmentCells[i],
        catchmentAreaKm2: catchmentAreaKm2[i],
        basinFloorAreaKm2: b.areaKm2,
        catchmentToFloorRatio: b.areaKm2 > 0 ? catchmentAreaKm2[i] / b.areaKm2 : null,
        drainTo,
        terminal,
    }));
    out.sharedRoots = sharedRoots;
    return out;
}

/**
 * Independently re-derive each basin's catchment area by summing the exported
 * `drainage_terminal` field, and check it against the per-basin areas. A
 * consumer can only trust the two products together if they agree, so the
 * agreement is measured and published rather than asserted in a comment.
 */
export function auditCatchmentConsistency(terminal, catchments, cellArea, opts = {}) {
    const tolerance = opts.tolerance ?? 1e-3;
    const byRoot = new Map(catchments.map(c => [c.sink, c]));
    const summed = new Map(catchments.map(c => [c.sink, 0]));
    for (let r = 0; r < terminal.length; r++) {
        const t = terminal[r];
        if (t >= 0 && summed.has(t)) summed.set(t, summed.get(t) + cellArea[r]);
    }
    let worst = 0, disagreeing = 0, totalSummed = 0, totalDeclared = 0;
    for (const [root, c] of byRoot) {
        const s = summed.get(root);
        totalSummed += s;
        totalDeclared += c.catchmentAreaKm2;
        const rel = c.catchmentAreaKm2 > 0 ? Math.abs(s - c.catchmentAreaKm2) / c.catchmentAreaKm2 : 0;
        if (rel > tolerance) disagreeing++;
        if (rel > worst) worst = rel;
    }
    return {
        basins: byRoot.size,
        disagreeing,
        worstRelativeError: worst,
        totalByFieldKm2: totalSummed,
        totalDeclaredKm2: totalDeclared,
        agrees: disagreeing === 0,
        // Publish the fraction WITH its denominator named. The totals above are
        // absolute km2 and denominator-free, which is safe but invites the
        // reader to divide by land_mask — the elevation-sign definition, which
        // drops the dry sub-sea-level basin floors and inflates this by ~3.5
        // points. A headline percentage is the number that gets quoted onward,
        // so it ships pre-computed against the right set.
        landAreaKm2: opts.landAreaKm2 ?? null,
        landDenominator: opts.landDenominator ?? null,
        fractionOfLand: opts.landAreaKm2 > 0 ? totalDeclared / opts.landAreaKm2 : null,
        denominatorWarning: 'Do NOT divide by the land_mask area. See manifest.landSeaMask: this '
            + 'export has two land definitions and they differ by the dry closed-basin floors, '
            + 'which are entirely endorheic and so bias this fraction upward when dropped.',
        note: 'Summing cell_area over drainage_terminal grouped by each basin\'s finalPreserved.sink '
            + 'must reproduce finalCatchment.areaKm2 exactly. Both are derived from the same routing; '
            + 'this is the published proof they were not derived from different sink definitions.',
    };
}

/**
 * Compare two catchment sets cell-for-cell — used to quantify how much erosion
 * moved a divide between the pre-conditioning and final terrain.
 */
export function compareCatchments(aCells, bCells, numRegions) {
    const inA = new Uint8Array(numRegions);
    for (const r of aCells) inA[r] = 1;
    let shared = 0;
    for (const r of bCells) if (inA[r]) shared++;
    const union = aCells.length + bCells.length - shared;
    return {
        aCount: aCells.length,
        bCount: bCells.length,
        sharedCount: shared,
        jaccard: union > 0 ? shared / union : 1,
        addedCount: bCells.length - shared,
        removedCount: aCells.length - shared,
    };
}

// ─────────────────────────────────────────────────────────────────────────
//  Remeasurement
// ─────────────────────────────────────────────────────────────────────────

/**
 * Re-measure a selected basin against the finished terrain.
 *
 * The headline number is spillDepthKm: how far below its spill point the sink
 * still sits. If protection worked, it stays close to the natural value; if the
 * carve passes got through, it collapses toward zero. This is the regression
 * signal for the whole feature.
 */
export function remeasureBasins(mesh, r_elevation, r_isOcean, selected, cellArea, reliefScale = 1) {
    const { numRegions, adjOffset, adjList } = mesh;
    const area = cellArea || uniformCellArea(numRegions);

    return selected.map(b => {
        const isMember = new Uint8Array(numRegions);
        for (const r of b.members) isMember[r] = 1;

        let sinkElev = Infinity, sink = b.sink;
        for (const r of b.members) {
            if (r_elevation[r] < sinkElev) { sinkElev = r_elevation[r]; sink = r; }
        }
        // Lowest saddle out of the member set on the current terrain.
        let spill = Infinity, spillFrom = -1, spillTo = -1;
        for (const r of b.members) {
            for (let i = adjOffset[r], iEnd = adjOffset[r + 1]; i < iEnd; i++) {
                const nb = adjList[i];
                if (isMember[nb]) continue;
                const saddle = Math.max(r_elevation[r], r_elevation[nb]);
                if (saddle < spill) { spill = saddle; spillFrom = r; spillTo = nb; }
            }
        }
        let volume = 0, floorArea = 0;
        if (spill < Infinity) {
            for (const r of b.members) {
                const d = spill - r_elevation[r];
                if (d > 0) { volume += d * area[r]; floorArea += area[r]; }
            }
        }
        const sinkH = heightKm(sinkElev, reliefScale);
        const spillH = spill === Infinity ? null : heightKm(spill, reliefScale);
        let volumeKm3 = 0;
        if (spill < Infinity) {
            for (const r of b.members) {
                const d = spillH - heightKm(r_elevation[r], reliefScale);
                if (d > 0) volumeKm3 += d * area[r];
            }
        }
        return {
            id: b.id,
            sink,
            sinkElevation: sinkElev,
            spillElevation: spill === Infinity ? null : spill,
            spillDepth: spill === Infinity ? null : spill - sinkElev,
            sinkElevationKm: sinkH,
            spillElevationKm: spillH,
            spillDepthKm: spillH === null ? null : spillH - sinkH,
            naturalSpillDepthKm: b.depthKm,
            // Retention compares like with like, in model units, so it is not
            // distorted by the nonlinear height curve.
            retainedFraction: b.depth > 0 && spill < Infinity
                ? (spill - sinkElev) / b.depth : null,
            volumeKm3,
            volumeModelUnits: volume,
            floodedAreaKm2: floorArea,
            spillFrom,
            spillTo,
        };
    });
}

// ─────────────────────────────────────────────────────────────────────────
//  Water levels (externally supplied — never decided here)
// ─────────────────────────────────────────────────────────────────────────

/**
 * Attach a lake surface elevation to preserved basins.
 *
 * Levels come from the caller — a downstream water balance, a hand-authored
 * worldbuilding decision, whatever. Orogen has no opinion. A basin with no
 * level supplied is a dry closed basin.
 *
 * Returns per-region surface classes and the validated level table.
 */
export function applyInlandWaterLevels(mesh, r_elevation, selected, levelsById = {}) {
    const { numRegions } = mesh;
    const waterLevelKm = new Float32Array(numRegions).fill(NaN);
    const inlandWaterBasin = new Int32Array(numRegions).fill(-1);
    const isInlandWater = new Uint8Array(numRegions);
    const warnings = [];
    const applied = [];

    for (const id of Object.keys(levelsById)) {
        if (!selected.some(b => b.id === id)) {
            warnings.push(`Water level given for "${id}", which is not a preserved basin — ignored.`);
        }
    }

    for (let i = 0; i < selected.length; i++) {
        const b = selected[i];
        const level = levelsById[b.id];
        if (level === undefined || level === null) { applied.push({ id: b.id, levelKm: null }); continue; }
        if (!Number.isFinite(level)) {
            warnings.push(`Water level for "${b.id}" is not a finite number — ignored.`);
            applied.push({ id: b.id, levelKm: null });
            continue;
        }
        if (level < b.sinkElevationKm) {
            warnings.push(
                `Water level ${level} km for "${b.id}" is below its floor ` +
                `(${b.sinkElevationKm.toFixed(4)} km) — the basin would be dry. Ignored.`);
            applied.push({ id: b.id, levelKm: null });
            continue;
        }
        if (level > b.spillElevationKm) {
            warnings.push(
                `Water level ${level} km for "${b.id}" is above its spill point ` +
                `(${b.spillElevationKm.toFixed(4)} km) — the basin would overflow and stop ` +
                `being endorheic. Clamped to the spill level.`);
        }
        const effective = Math.min(level, b.spillElevationKm);
        let area = 0, volume = 0;
        for (const r of b.members) {
            if (r_elevation[r] <= effective) {
                isInlandWater[r] = 1;
                inlandWaterBasin[r] = i;
                waterLevelKm[r] = effective;
                area += 1;
                volume += effective - r_elevation[r];
            }
        }
        applied.push({ id: b.id, levelKm: effective, floodedCells: area, meanDepthKm: area > 0 ? volume / area : 0 });
    }

    return { isInlandWater, inlandWaterBasin, waterLevelKm, applied, warnings };
}

/**
 * Surface class per region: 0 ocean, 1 land, 2 inland water.
 *
 * Takes the OPEN-ocean mask, not an elevation<=0 test. A closed basin floor
 * below sea level is not the sea — classing it as ocean is exactly the mistake
 * this module exists to prevent, and it silently swallowed the inland-water
 * cells of every below-sea-level preserved basin. A dry closed basin below sea
 * level is land; it becomes inland water only when a level is supplied.
 */
export function surfaceClasses(r_elevation, isOpenOcean, isInlandWater) {
    const n = r_elevation.length;
    const out = new Uint8Array(n);
    for (let r = 0; r < n; r++) {
        if (isInlandWater && isInlandWater[r]) out[r] = 2;
        else out[r] = isOpenOcean[r] ? 0 : 1;
    }
    return out;
}

export const BASIN_SURFACE_CLASSES = { 0: 'ocean', 1: 'land', 2: 'inland_water' };

/**
 * Surface classes whose elevation converts through the LAND branch of the height
 * curve — subaerial ground and the floor beneath an inland lake, both of which
 * are continental crust, against open seabed which is not.
 *
 * Stated as an explicit set rather than `!== ocean`. A set defined as the
 * complement of another is not wrong when written and does not go stale in any
 * way a schema or hash check can see — it breaks silently the moment someone
 * adds a member to the vocabulary it negates, and the new member joins the
 * complement by default. Downstream lost a whole rock class that way: their
 * "vegetated" mode meant `rock != evaporite`, which was correct until evaporite
 * was split and playa fill silently inherited a canopy. Say the set.
 */
export const LAND_HEIGHT_BRANCH_CLASSES = Object.freeze([1, 2]);

/** True when this surface class uses the land branch of the height curve. */
export function usesLandHeightBranch(surfaceClass) {
    return LAND_HEIGHT_BRANCH_CLASSES.includes(surfaceClass);
}

/** Resolution context, so a catalogue can be read without guessing the mesh. */
export function basinResolutionContext(numRegions, radiusKm, opts = {}) {
    // Report the thresholds ACTUALLY used, not the config defaults. Reading the
    // constants here meant a run with overridden thresholds published the wrong
    // numbers, which is worse than publishing none — a downstream consumer would
    // reason about a selection that never happened.
    const minDepthKm = opts.minDepthKm ?? BASIN_MIN_DEPTH_KM;
    const minAreaKm2 = opts.minAreaKm2 ?? BASIN_MIN_AREA_KM2;
    const minCells = opts.minCells ?? BASIN_MIN_CELLS;
    const cell = avgCellAreaKm2(numRegions, radiusKm);
    return {
        numRegions,
        radiusKm: radiusKm ?? null,
        avgEdgeKm: avgEdgeKm(numRegions, radiusKm),
        avgCellAreaKm2: cell,
        minDepthKm,
        minAreaKm2,
        minCells,
        // The currency the depth floor is compared in, stated rather than left
        // to be inferred: `minDepthKm` above is enforced against the
        // depression's physical depth, so the declared number is the drop it
        // demands everywhere on the height curve. Published because a consumer
        // reading a criteria block cannot otherwise tell a physical floor from
        // a model-unit one wearing a physical name, and the two select very
        // different catalogues. `selectBasins` has the argument.
        minDepthComparedIn: 'km',
        minDepthIsPhysical: true,
        // Which floor binds, ON THE AVERAGE CELL. `cell` is the mean dual area
        // over the whole sphere, so this is a statement about the typical
        // depression and not a guarantee about any particular one: local cell
        // area varies with the mesh jitter, and where cells run large the cell
        // floor can still reject a depression above effectiveMinAreaKm2. Ask
        // `hydrography/scripts/catalogue_floor.py` what it did on a given
        // catalogue rather than inferring it from this field.
        effectiveMinAreaKm2: Math.max(minAreaKm2, minCells * cell),
        bindingFloor: minCells * cell > minAreaKm2 ? 'minCells' : 'minAreaKm2',
        bindingFloorBasis: 'mean cell area over the sphere',
        // And neither: the two floors above compare area against area, while
        // the depth floor compares a length, so `bindingFloor` never names it
        // however much of the catalogue it decides. Ask
        // `hydrography/scripts/catalogue_floor.py`, which relaxes each floor and
        // re-runs the whole selection, which of the three actually decided a
        // given build.
        bindingFloorExcludes: 'minDepthKm',
    };
}

// ─────────────────────────────────────────────────────────────────────────
//  Rivers
// ─────────────────────────────────────────────────────────────────────────

/**
 * Upstream drainage area and river mouths on the finished terrain.
 *
 * erodeComposite builds a flow-accumulation array on every hydraulic iteration
 * and throws it away. This recomputes it once against the final surface, where
 * it means something stable, and in real units: upstream contributing area in
 * km², not a cell count, so it compares across resolutions and planets.
 *
 * River mouths are where a land cell discharges into the open ocean. Their
 * accumulated area is the catchment feeding that point, which is what a
 * downstream model wants for freshwater flux — the single largest control on
 * coastal salinity, and the thing that tells a worldbuilder where the great
 * rivers reach the sea.
 */
export function computeRivers(mesh, r_elevation, r_isOcean, cellArea, opts = {}) {
    const { numRegions, adjOffset, adjList } = mesh;
    const area = cellArea || uniformCellArea(numRegions, opts.radiusKm);
    const terminalRoot = opts.terminalRoot || null;
    const { isOpenOcean } = findOpenOcean(mesh, r_isOcean);

    // ONE TREE, THREE PRODUCTS. drain_to, flow_accum_km2 and drainage_terminal
    // are three views of the same routing, so they come from the same call. They
    // used to be built separately, which produced two failures at once:
    //
    //   - This function did raw steepest descent, so every cell draining into a
    //     post-erosion pit accumulated into that pit and stopped.
    //   - It was handed the basin MEMBER mask as terminalRoot, not the sink set,
    //     so a river reaching a preserved basin terminated at the rim it crossed
    //     and the sink itself showed only its own cell area — no lake inflow
    //     anywhere on the planet.
    //
    // computeFinalTerrainDrainage routes on the filled surface via the flood's
    // visitation tree, which is exactly steepest descent on that surface (the
    // heap pops by fill value, so the neighbour that first reaches a cell is its
    // lowest-fill neighbour) and additionally crosses filled flats, where
    // steepest descent has no gradient at all.
    const { drainTo } = computeFinalTerrainDrainage(
        mesh, r_elevation, r_isOcean, terminalRoot, { raw: !!opts.raw });

    // Accumulate downhill over that tree. Topological order by in-degree rather
    // than by elevation: on a filled flat, many cells share an elevation and a
    // sort gives no guarantee that a cell is complete before it is passed on.
    const accum = new Float32Array(numRegions);
    for (let r = 0; r < numRegions; r++) accum[r] = area[r];
    const indeg = new Int32Array(numRegions);
    for (let r = 0; r < numRegions; r++) if (drainTo[r] >= 0) indeg[drainTo[r]]++;
    const queue = new Int32Array(numRegions);
    let head = 0, tail = 0;
    for (let r = 0; r < numRegions; r++) if (indeg[r] === 0) queue[tail++] = r;
    while (head < tail) {
        const r = queue[head++];
        const t = drainTo[r];
        if (t < 0) continue;
        accum[t] += accum[r];
        if (--indeg[t] === 0) queue[tail++] = t;
    }

    // Mouths: land cells that discharge straight into the open ocean.
    const mouths = [];
    for (let r = 0; r < numRegions; r++) {
        if (isOpenOcean[r] || r_elevation[r] <= 0) continue;
        const t = drainTo[r];
        if (t >= 0 && isOpenOcean[t]) {
            mouths.push({ region: r, dischargeAreaKm2: accum[r] });
        }
    }
    mouths.sort((a, b) => b.dischargeAreaKm2 - a.dischargeAreaKm2 || a.region - b.region);

    return { flowAccumKm2: accum, drainTo, mouths, isOpenOcean };
}

// ─────────────────────────────────────────────────────────────────────────
//  Ocean basins and sills
// ─────────────────────────────────────────────────────────────────────────

/**
 * Marginal seas, and the sill depth separating each from the world ocean.
 *
 * Exactly the closed-basin problem, run below sea level instead of above it: a
 * marginal sea is a depression in the bathymetry whose rim is the sill. The
 * sill depth is what decides the character of the basin — a deep sill exchanges
 * freely with the world ocean, a shallow one restricts it into a Mediterranean
 * (evaporitic, salty) or a Black Sea (stratified, anoxic below the pycnocline).
 * That distinction matters both to an ocean model and to worldbuilding, and it
 * is not recoverable from bathymetry alone without this computation.
 */
export function computeOceanBasins(mesh, r_elevation, r_isOcean, cellArea, opts = {}) {
    const { numRegions, adjOffset, adjList } = mesh;
    const area = cellArea || uniformCellArea(numRegions, opts.radiusKm);
    const { isOpenOcean } = findOpenOcean(mesh, r_isOcean);

    // Connected components of submerged cells that are NOT the world ocean.
    const label = new Int32Array(numRegions).fill(-1);
    const basins = [];
    for (let seed = 0; seed < numRegions; seed++) {
        if (label[seed] >= 0 || isOpenOcean[seed] || r_elevation[seed] > 0) continue;
        const id = basins.length;
        const members = [];
        const stack = [seed];
        label[seed] = id;
        while (stack.length > 0) {
            const r = stack.pop();
            members.push(r);
            for (let i = adjOffset[r], iEnd = adjOffset[r + 1]; i < iEnd; i++) {
                const nb = adjList[i];
                if (label[nb] >= 0 || isOpenOcean[nb] || r_elevation[nb] > 0) continue;
                label[nb] = id;
                stack.push(nb);
            }
        }
        members.sort((a, b) => a - b);

        // Sill: the lowest point on the rim, i.e. the shallowest barrier that
        // must be crossed to leave. Its elevation IS the sill depth (negative
        // below sea level); a sill above sea level means fully enclosed.
        let sillElev = Infinity, sillFrom = -1, sillTo = -1, touchesOpenOcean = false;
        let deepest = Infinity, deepestCell = members[0], areaKm2 = 0, volumeKm3 = 0;
        for (const r of members) {
            areaKm2 += area[r];
            if (r_elevation[r] < deepest) { deepest = r_elevation[r]; deepestCell = r; }
            for (let i = adjOffset[r], iEnd = adjOffset[r + 1]; i < iEnd; i++) {
                const nb = adjList[i];
                if (label[nb] === id) continue;
                if (isOpenOcean[nb]) touchesOpenOcean = true;
                const barrier = Math.max(r_elevation[r], r_elevation[nb]);
                if (barrier < sillElev) { sillElev = barrier; sillFrom = r; sillTo = nb; }
            }
        }
        for (const r of members) {
            const d = Math.min(0, sillElev) - r_elevation[r];
            if (d > 0) volumeKm3 += d * area[r];
        }

        basins.push({
            index: id,
            id: basinId(deepestCell, members),
            deepestCell,
            maxDepthKm: -deepest,
            sillDepthKm: sillElev === Infinity ? null : -Math.min(0, sillElev),
            sillElevationKm: sillElev === Infinity ? null : sillElev,
            sillFrom, sillTo,
            enclosed: sillElev > 0,
            adjacentToOpenOcean: touchesOpenOcean,
            areaKm2,
            volumeBelowSillKm3: volumeKm3,
            cellCount: members.length,
            members,
        });
    }
    basins.sort((a, b) => b.areaKm2 - a.areaKm2 || a.deepestCell - b.deepestCell);
    return { basins, label, isOpenOcean };
}
