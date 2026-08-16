#!/usr/bin/env node
/**
 * Headless equirectangular map-image export.
 *
 * The browser's exportMap() renders the dual mesh on the GPU. This is the CPU
 * equivalent: it rasterizes the same triangles, with the same colour functions
 * and the same antimeridian handling, so the output matches what the Export Map
 * button produces — without a browser, at any size, and scriptable.
 *
 * Requires: npm i delaunator
 *
 * Usage:
 *   node tools/export-maps.mjs --code 09oa7lj8kek17v5639gvhe --width 16384
 *   node tools/export-maps.mjs --seed 12345 --regions 250000 --width 4096 --types color,landmask
 *
 * Run with --help for the full option list.
 */

import fs from 'node:fs';
import path from 'node:path';
import zlib from 'node:zlib';
import Delaunator from 'delaunator';

import { setDelaunator } from '../js/sphere-mesh.js';
import { decodePlanetCode } from '../js/planet-code.js';
import { gravityFromMassRadius, makePlanet, planetSummary } from '../js/planet-params.js';
import { basinAreaFromSlider, BASIN_AREA_LADDER_KM2 } from '../js/terrain-config.js';
import { encodePlanetCode } from '../js/planet-code.js';

/** Nearest ladder index for a basin floor area, for re-encoding a planet code. */
function sliderFromBasinArea(km2) {
    if (!(km2 > 0)) return 0;
    let best = 0, bestD = Infinity;
    for (let i = 1; i < BASIN_AREA_LADDER_KM2.length; i++) {
        const d = Math.abs(BASIN_AREA_LADDER_KM2[i] - km2);
        if (d < bestD) { bestD = d; best = i; }
    }
    return best;
}
import { elevToHeightKm, scaledHeightKm, elevationToColor, biomeColor } from '../js/color-map.js';
import { KOPPEN_CLASSES } from '../js/koppen.js';
import { usesLandHeightBranch, parseBasinList } from '../js/basins.js';
import { hashTypedArray, hashJson } from '../js/sha256.js';

setDelaunator(Delaunator);

// ── Colour functions, verbatim from js/planet-mesh.js ────────────────────

// Physical height for this planet: land carries the 1/g relief scaling, ocean
// does not (g cancels in the isostatic balance of the water column).
// Physical height for this planet. `isLand` must come from surface_class, not
// the elevation sign: a dry closed-basin floor is land, and the bathymetric
// branch rendered ~-500 m floors at -4,300 to -4,900 m, clamping a quarter of a
// million pixels at the encoding floor of the full-range heightmap.
const heightKm = (e, reliefScale, isLand) => scaledHeightKm(e, reliefScale, isLand);
const heightmapColor = (e, reliefScale, isLand) => {
    const t = Math.max(0, Math.min(1, (heightKm(e, reliefScale, isLand) + 5) / 11));  // -5 km → 0, 6 km → 1
    return [t, t, t];
};
const landHeightmapColor = (e, reliefScale, isLand) => {
    if (e <= 0) return [0, 0, 0];
    const t = Math.max(0, Math.min(1, heightKm(e, reliefScale, isLand) / 6));
    return [t, t, t];
};
const landMaskColor = (e) => (e > 0 ? [1, 1, 1] : [0, 0, 0]);
const koppenColor = (id) => KOPPEN_CLASSES[id] ? KOPPEN_CLASSES[id].color : [0, 0, 0];

const TYPES = {
    color:          { file: 'colormap',        bw: false, gray16: false },
    biome:          { file: 'satellite',       bw: false, gray16: false, climate: true },
    koppen:         { file: 'climate',         bw: false, gray16: false, climate: true },
    heightmap:      { file: 'heightmap',       bw: true,  gray16: true  },
    landheightmap:  { file: 'land-heightmap',  bw: true,  gray16: true  },
    landmask:       { file: 'landmask',        bw: true,  gray16: false },
    surfacemask:    { file: 'surfacemask',     bw: true,  gray16: false },
};

// ── PNG encoding ─────────────────────────────────────────────────────────

const CRC_TABLE = (() => {
    const t = new Uint32Array(256);
    for (let n = 0; n < 256; n++) {
        let c = n;
        for (let k = 0; k < 8; k++) c = (c & 1) ? (0xEDB88320 ^ (c >>> 1)) : (c >>> 1);
        t[n] = c >>> 0;
    }
    return t;
})();

function crc32(buf) {
    let c = 0xFFFFFFFF;
    for (let i = 0; i < buf.length; i++) c = CRC_TABLE[(c ^ buf[i]) & 0xFF] ^ (c >>> 8);
    return (c ^ 0xFFFFFFFF) >>> 0;
}

function chunk(type, data) {
    const out = Buffer.alloc(12 + data.length);
    out.writeUInt32BE(data.length, 0);
    out.write(type, 4, 'ascii');
    data.copy(out, 8);
    out.writeUInt32BE(crc32(out.subarray(4, 8 + data.length)), 8 + data.length);
    return out;
}

/**
 * Write a PNG, streaming scanlines through deflate so the compressed image
 * never has to exist in memory alongside the raw one. A 16384×8192 RGBA raster
 * is 537 MB raw; holding a second copy of it would be gratuitous.
 *
 * colorType 6 = RGBA/8, 0 = grayscale (8 or 16 bit). Filter type 0 (None) on
 * every scanline — these images are mostly smooth gradients or flat regions and
 * the extra CPU of adaptive filtering is not worth it at this size.
 */
function writePNG(file, width, height, colorType, bitDepth, getRow) {
    const channels = colorType === 6 ? 4 : 1;
    const bytesPerPixel = channels * (bitDepth / 8);
    const stride = width * bytesPerPixel;

    const fd = fs.openSync(file, 'w');
    fs.writeSync(fd, Buffer.from([0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A]));

    const ihdr = Buffer.alloc(13);
    ihdr.writeUInt32BE(width, 0);
    ihdr.writeUInt32BE(height, 4);
    ihdr[8] = bitDepth;
    ihdr[9] = colorType;
    ihdr[10] = 0; ihdr[11] = 0; ihdr[12] = 0;   // deflate / adaptive filtering / no interlace
    fs.writeSync(fd, chunk('IHDR', ihdr));

    const deflate = zlib.deflateSync(
        (function () {
            // Assemble filtered scanlines once; deflateSync over the whole thing
            // is markedly faster than streaming for this shape of data.
            const buf = Buffer.alloc((stride + 1) * height);
            let off = 0;
            for (let y = 0; y < height; y++) {
                buf[off++] = 0;                      // filter: None
                getRow(y, buf, off);
                off += stride;
            }
            return buf;
        })(),
        { level: 6 }
    );
    fs.writeSync(fd, chunk('IDAT', deflate));
    fs.writeSync(fd, chunk('IEND', Buffer.alloc(0)));
    fs.closeSync(fd);
    return fs.statSync(file).size;
}

// ── Rasterizer ───────────────────────────────────────────────────────────

/**
 * Rasterize the dual mesh into an equirectangular raster.
 *
 * Mirrors buildMapMesh/exportMap: one triangle per half-edge, spanning
 * (inner triangle centre, outer triangle centre, region centre). Vertex colours
 * are interpolated barycentrically, which reproduces the GPU's Gouraud shading
 * for the smooth heightmap types; for flat types all three vertices carry the
 * same colour so interpolation is a no-op.
 *
 * Antimeridian: a triangle whose longitude span exceeds π is emitted twice,
 * shifted by ±2π, exactly as the browser does.
 */
function rasterize(mesh, r_xyz, t_xyz, width, height, vertexColor, channels, background, gray16, onProgress) {
    const maxV = gray16 ? 65535 : 255;
    const img = gray16 ? new Uint16Array(width * height * channels)
                       : new Uint8Array(width * height * channels);
    for (let i = 0; i < img.length; i += channels) {
        for (let c = 0; c < channels; c++) {
            img[i + c] = Math.max(0, Math.min(maxV, Math.round(background[c] * maxV)));
        }
    }

    const PI = Math.PI;
    const clamp = (v, lo, hi) => (v < lo ? lo : v > hi ? hi : v);
    // lon (rad) → pixel x, lat (rad) → pixel y. Cell-centre registered, so a
    // point at lon = -π lands on the left edge of column 0.
    const px = (lon) => (lon / (2 * PI) + 0.5) * width - 0.5;
    const py = (lat) => (0.5 - lat / PI) * height - 0.5;

    const cx = new Float64Array(3), cy = new Float64Array(3);
    const col = [new Float64Array(channels), new Float64Array(channels), new Float64Array(channels)];

    const drawTriangle = () => {
        let minX = Math.floor(Math.min(cx[0], cx[1], cx[2]));
        let maxX = Math.ceil(Math.max(cx[0], cx[1], cx[2]));
        let minY = Math.floor(Math.min(cy[0], cy[1], cy[2]));
        let maxY = Math.ceil(Math.max(cy[0], cy[1], cy[2]));
        if (maxX < 0 || minX >= width || maxY < 0 || minY >= height) return;
        minX = clamp(minX, 0, width - 1); maxX = clamp(maxX, 0, width - 1);
        minY = clamp(minY, 0, height - 1); maxY = clamp(maxY, 0, height - 1);

        const x0 = cx[0], y0 = cy[0], x1 = cx[1], y1 = cy[1], x2 = cx[2], y2 = cy[2];
        const denom = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2);
        if (Math.abs(denom) < 1e-12) return;   // degenerate
        const inv = 1 / denom;

        for (let y = minY; y <= maxY; y++) {
            const fy = y - y2;
            for (let x = minX; x <= maxX; x++) {
                const fx = x - x2;
                const l0 = ((y1 - y2) * fx + (x2 - x1) * fy) * inv;
                if (l0 < -1e-9) continue;
                const l1 = ((y2 - y0) * fx + (x0 - x2) * fy) * inv;
                if (l1 < -1e-9) continue;
                const l2 = 1 - l0 - l1;
                if (l2 < -1e-9) continue;
                const o = (y * width + x) * channels;
                for (let c = 0; c < channels; c++) {
                    const v = (l0 * col[0][c] + l1 * col[1][c] + l2 * col[2][c]) * maxV;
                    img[o + c] = v < 0 ? 0 : v > maxV ? maxV : Math.round(v);
                }
            }
        }
    };

    const numSides = mesh.numSides;
    const reportEvery = Math.max(1, Math.floor(numSides / 20));

    for (let s = 0; s < numSides; s++) {
        if (onProgress && s % reportEvery === 0) onProgress(s / numSides);

        const it = mesh.s_inner_t(s), ot = mesh.s_outer_t(s), br = mesh.s_begin_r(s);
        vertexColor(s, it, ot, br, col);

        const ax = t_xyz[3 * it], ay = t_xyz[3 * it + 1], az = t_xyz[3 * it + 2];
        const bx = t_xyz[3 * ot], by = t_xyz[3 * ot + 1], bz = t_xyz[3 * ot + 2];
        const gx = r_xyz[3 * br], gy = r_xyz[3 * br + 1], gz = r_xyz[3 * br + 2];

        let lon0 = Math.atan2(ax, az), lat0 = Math.asin(clamp(ay, -1, 1));
        let lon1 = Math.atan2(bx, bz), lat1 = Math.asin(clamp(by, -1, 1));
        let lon2 = Math.atan2(gx, gz), lat2 = Math.asin(clamp(gy, -1, 1));

        const wraps = Math.max(lon0, lon1, lon2) - Math.min(lon0, lon1, lon2) > PI;
        if (wraps) {
            if (lon0 < 0) lon0 += 2 * PI;
            if (lon1 < 0) lon1 += 2 * PI;
            if (lon2 < 0) lon2 += 2 * PI;
        }

        cy[0] = py(lat0); cy[1] = py(lat1); cy[2] = py(lat2);
        cx[0] = px(lon0); cx[1] = px(lon1); cx[2] = px(lon2);
        drawTriangle();

        if (wraps) {
            cx[0] = px(lon0 - 2 * PI); cx[1] = px(lon1 - 2 * PI); cx[2] = px(lon2 - 2 * PI);
            drawTriangle();
        }
    }
    if (onProgress) onProgress(1);
    return img;
}

// ── CLI ──────────────────────────────────────────────────────────────────

const HELP = `
Headless equirectangular map-image export.

  --code STR       planet code; supplies seed and every slider
  --seed N         seed (ignored when --code is given)
  --regions N      mesh region count (ignored when --code is given)
  --width N        image width in pixels; height is width/2 (default: 4096)
  --types A,B      any of: ${Object.keys(TYPES).join(', ')}
                   landmask is the elevation-sign convention and floods dry closed
                   basins below sea level; surfacemask is the real land/sea mask
                   (white land, grey inland water, black ocean). Prefer surfacemask.
                   (default: color,heightmap,landheightmap,landmask)
  --out DIR        output directory (default: out/maps-<code|seed>)
  --name STR       stem used in the filenames (default: the code, else the seed)
  --no-basins      disable endorheic basin preservation (upstream drainage)
  --no-lithology   uniform erodibility (upstream erosion, no differential relief)
  --basins N       endorheic-basin slider index 0-8 (0 = off); overrides the code
  --basin-min-area KM2   basin floor-area threshold, overriding the slider ladder
  --basin-min-cells N    basin resolution floor in mesh cells
  --basin-min-depth KM   basin depth-below-spill threshold
  --preserve-basins FILE explicit drainage hypothesis, same format and meaning as
                         export-planet.mjs: the preserved set becomes exactly the
                         basins listed, ignoring the size floors. Needed to render
                         a map of a data export built from a carve verdict — the
                         planet code cannot carry one.
  --carve-basins FILE    same format, subtractive: threshold selection minus these.
  --rock-contrast F  lithology strength 0-1; overrides the code
  --radius KM      planetary radius (default 6371)
  --gravity MS2    surface gravity (default 9.80665)
  --mass-radius M,R  gravity from Earth masses and Earth radii instead
  --planet-name S  label for the planet
  --quiet
  --help

Naming matches the browser's Export Map button: orogen-<type>-<name>.png
`;

function parseArgs(argv) {
    const a = {
        code: null, seed: null, regions: null, width: 4096,
        types: ['color', 'heightmap', 'landheightmap', 'landmask', 'surfacemask'],
        out: null, name: null, quiet: false, preserveBasins: true, lithology: true,
        planet: {}, basins: null, rockContrast: null,
        basinMinAreaKm2: null, basinMinCells: null, basinMinDepthKm: null,
        preserveBasinList: null, carveBasinList: null, basinListFile: null,
    };
    const num = (v, f) => {
        const n = Number(v);
        if (!Number.isFinite(n)) throw new Error(`${f} expects a number, got "${v}"`);
        return n;
    };
    for (let i = 2; i < argv.length; i++) {
        const f = argv[i];
        switch (f) {
            case '--help': case '-h': console.log(HELP); process.exit(0); break;
            case '--code': a.code = argv[++i]; break;
            case '--seed': a.seed = num(argv[++i], f); break;
            case '--regions': case '--n': a.regions = num(argv[++i], f); break;
            case '--width': a.width = num(argv[++i], f); break;
            case '--types': a.types = argv[++i].split(',').map(s => s.trim()).filter(Boolean); break;
            case '--out': a.out = argv[++i]; break;
            case '--name': a.name = argv[++i]; break;
            case '--no-basins': a.preserveBasins = false; break;
            case '--no-lithology': a.lithology = false; break;
            case '--basins': a.basins = num(argv[++i], f); break;
            case '--basin-min-area': a.basinMinAreaKm2 = num(argv[++i], f); break;
            case '--basin-min-cells': a.basinMinCells = num(argv[++i], f); break;
            case '--basin-min-depth': a.basinMinDepthKm = num(argv[++i], f); break;
            case '--preserve-basins': a.basinListFile = argv[++i]; a.preserveBasinList = parseBasinList(fs.readFileSync(a.basinListFile, 'utf8')); break;
            case '--carve-basins': a.basinListFile = argv[++i]; a.carveBasinList = parseBasinList(fs.readFileSync(a.basinListFile, 'utf8')); break;
            case '--rock-contrast': a.rockContrast = num(argv[++i], f); break;
            case '--radius': a.planet.radiusKm = num(argv[++i], f); break;
            case '--gravity': a.planet.gravityMS2 = num(argv[++i], f); break;
            case '--planet-name': a.planet.name = argv[++i]; break;
            case '--mass-radius': {
                const parts = String(argv[++i]).split(',');
                if (parts.length !== 2) throw new Error('--mass-radius expects M,R in Earth units');
                a.planet.gravityMS2 = gravityFromMassRadius(Number(parts[0]), Number(parts[1]));
                break;
            }
            case '--quiet': a.quiet = true; break;
            default: throw new Error(`Unknown option: ${f} (try --help)`);
        }
    }
    const bad = a.types.filter(t => !TYPES[t]);
    if (bad.length) throw new Error(`Unknown type(s): ${bad.join(', ')}. Valid: ${Object.keys(TYPES).join(', ')}`);
    if (!a.code && a.seed === null) throw new Error('Give either --code or --seed');
    if (a.width < 2 || (a.width & 1)) throw new Error('--width must be an even number >= 2');
    return a;
}

function humanBytes(n) {
    const u = ['B', 'KB', 'MB', 'GB'];
    let i = 0;
    while (n >= 1024 && i < u.length - 1) { n /= 1024; i++; }
    return `${n.toFixed(i === 0 ? 0 : 1)} ${u[i]}`;
}

async function main() {
    const args = parseArgs(process.argv);
    const log = args.quiet ? () => {} : (...m) => console.log(...m);
    const stdout = console.log.bind(console);

    let params;
    if (args.code) {
        const d = decodePlanetCode(args.code);
        if (!d) throw new Error(`Could not decode planet code "${args.code}"`);
        params = {
            N: d.N, P: d.P, jitter: d.jitter, nMag: d.roughness,
            numContinents: d.numContinents, continentSizeVariety: d.continentSizeVariety,
            landCoverage: d.landCoverage, terrainWarp: d.terrainWarp, smoothing: d.smoothing,
            hydraulicErosion: d.hydraulicErosion, thermalErosion: d.thermalErosion,
            ridgeSharpening: d.ridgeSharpening, glacialErosion: d.glacialErosion,
            temperatureOffset: d.temperatureOffset, precipitationOffset: d.precipitationOffset,
            seed: d.seed, toggledIndices: d.toggledIndices || [],
        };
        // Parity with export-planet.mjs: a code carries basin and rock-contrast
        // settings and they must be honoured, not silently dropped. Explicit
        // flags still win — handled below.
        const basinKm2 = basinAreaFromSlider(d.basinSlider ?? 0);
        params.preserveBasins = basinKm2 > 0;
        if (basinKm2 > 0) params.basinMinAreaKm2 = basinKm2;
        const strength = d.lithologyStrength ?? 0;
        params.lithology = strength > 0;
        params.lithologyStrength = strength;
        if (args.regions !== null) params.N = args.regions;
    } else {
        params = {
            N: args.regions ?? 100000, P: 24, jitter: 0.75, nMag: 1.0,
            numContinents: 5, continentSizeVariety: 0, landCoverage: 0.3,
            terrainWarp: 0, smoothing: 0, hydraulicErosion: 0.4, thermalErosion: 0.3,
            ridgeSharpening: 0.3, glacialErosion: 0.3,
            temperatureOffset: 0, precipitationOffset: 0,
            seed: args.seed, toggledIndices: [],
        };
    }

    // Name the files after the planet that was actually generated. Naming them
    // after the input code is a trap: any override (--regions above all, since
    // region count is a generation parameter) makes the code describe a
    // different terrain than the files contain.
    let effectiveCode = null;
    try {
        effectiveCode = encodePlanetCode(
            params.seed, params.N, params.jitter, params.P, params.numContinents, params.nMag,
            params.terrainWarp, params.smoothing, params.glacialErosion, params.hydraulicErosion,
            params.thermalErosion, params.ridgeSharpening, 0.75, params.continentSizeVariety,
            params.temperatureOffset, params.precipitationOffset, params.landCoverage,
            sliderFromBasinArea(params.preserveBasins ? (params.basinMinAreaKm2 ?? 1000) : 0),
            params.lithology ? (params.lithologyStrength ?? 1) : 0,
            params.toggledIndices);
    } catch { /* out-of-range slider: fall back to the seed */ }

    if (args.code && effectiveCode && effectiveCode !== args.code && !args.name) {
        console.warn(`  note: overrides changed the planet; naming output after the effective code`);
        console.warn(`        given ${args.code} -> effective ${effectiveCode}`);
    }
    if ((args.preserveBasinList || args.carveBasinList) && !args.name) {
        console.warn(`  note: a basin list changes the terrain and no planet code can express it;`);
        console.warn(`        the filenames say ${effectiveCode ?? params.seed}, which is this planet WITHOUT the list. Pass --name.`);
    }
    const name = args.name ?? effectiveCode ?? String(params.seed);
    const outDir = args.out ?? path.join('out', `maps-${name}`);
    const needClimate = args.types.some(t => TYPES[t].climate);
    params.skipClimate = !needClimate;
    params.preserveBasins = args.preserveBasins;
    params.lithology = args.lithology;
    params.planet = args.planet;
    if (args.basins !== null) {
        const km2 = basinAreaFromSlider(args.basins);
        params.preserveBasins = km2 > 0;
        if (km2 > 0) params.basinMinAreaKm2 = km2;
    }
    // Any explicit threshold means basins are wanted, whatever the slider says.
    if (args.basinMinAreaKm2 !== null) { params.preserveBasins = true; params.basinMinAreaKm2 = args.basinMinAreaKm2; }
    if (args.basinMinCells !== null) { params.preserveBasins = true; params.basinMinCells = args.basinMinCells; }
    if (args.basinMinDepthKm !== null) { params.preserveBasins = true; params.basinMinDepthKm = args.basinMinDepthKm; }
    // A drainage hypothesis from downstream, as export-planet.mjs takes it. A
    // planet code cannot carry one, so rendering this planet without the list
    // that its data export was built from silently draws a different terrain.
    if (args.preserveBasinList || args.carveBasinList) params.preserveBasins = true;
    params.preserveBasinList = args.preserveBasinList;
    params.carveBasinList = args.carveBasinList;
    if (args.rockContrast !== null) {
        params.lithology = args.rockContrast > 0;
        params.lithologyStrength = args.rockContrast;
    }
    const pl = makePlanet(args.planet);
    if (pl.radiusKm !== 6371 || pl.gravityMS2 !== 9.80665) {
        log(`  ${pl.name}: R=${pl.radiusKm} km (${pl.radiusEarth.toFixed(2)} Earth), `
          + `g=${pl.gravityMS2.toFixed(3)} m/s² (${pl.gravityEarth.toFixed(2)} Earth), `
          + `relief x${pl.reliefScale.toFixed(3)}`);
    }

    const width = args.width, height = width / 2;
    log(`Planet ${name}: seed ${params.seed}, ${params.N.toLocaleString()} regions, ${params.P} plates`);
    log(`Rendering ${args.types.join(', ')} at ${width}×${height}${needClimate ? ' (climate on)' : ''}`);

    if (args.quiet) console.log = () => {};
    const { runGeneratePipeline } = await import('../js/pipeline.js');

    const t0 = performance.now();
    let lastPct = -1;
    const ctx = runGeneratePipeline(params, (pct, label) => {
        if (pct !== lastPct) { log(`  [${String(pct).padStart(3)}%] ${label}`); lastPct = pct; }
    });
    log(`Generated in ${((performance.now() - t0) / 1000).toFixed(1)}s`);
    if (ctx.basins) log(`Basins: ${ctx.basins.catalogue.length} depressions, ${ctx.basins.selected.length} preserved`);

    const { mesh, r_xyz, t_xyz, r_elevation } = ctx;
    const koppen = ctx.debugLayers ? ctx.debugLayers.koppen : null;
    const surfaceClass = ctx.basins ? ctx.basins.r_surfaceClass : null;

    // Triangle-centre elevations, so the smooth heightmap types interpolate the
    // same values the GPU path does.
    const t_elev = new Float32Array(mesh.numTriangles);
    for (let t = 0; t < mesh.numTriangles; t++) {
        const s0 = 3 * t;
        t_elev[t] = (r_elevation[mesh.s_begin_r(s0)] + r_elevation[mesh.s_begin_r(s0 + 1)]
                   + r_elevation[mesh.s_begin_r(s0 + 2)]) / 3;
    }

    fs.mkdirSync(outDir, { recursive: true });
    const written = [];

    for (const type of args.types) {
        const spec = TYPES[type];
        if (spec.climate && !koppen) {
            console.warn(`  skipping ${type}: needs climate, which was not computed`);
            continue;
        }
        const channels = spec.gray16 ? 1 : 3;
        const background = spec.bw ? [0, 0, 0] : [0x1a / 255, 0x1a / 255, 0x2e / 255];

        let vertexColor;
        if (type === 'heightmap' || type === 'landheightmap') {
            // Smooth: triangle-centre vertices carry averaged elevation.
            const fn = type === 'landheightmap' ? landHeightmapColor : heightmapColor;
            const rs = ctx.planet ? ctx.planet.reliefScale : 1;
            // Triangle-centre vertices take the land flag of the region they
            // belong to; a triangle straddling a basin rim interpolates between
            // the two, which is what Gouraud shading is for.
            const isLandOf = (r) => (surfaceClass ? usesLandHeightBranch(surfaceClass[r]) : r_elevation[r] > 0);
            vertexColor = (s, it, ot, br, col) => {
                const L = isLandOf(br);
                col[0][0] = fn(t_elev[it], rs, L)[0];
                col[1][0] = fn(t_elev[ot], rs, L)[0];
                col[2][0] = fn(r_elevation[br], rs, L)[0];
            };
        } else {
            const fn = type === 'surfacemask' ? (e, r) => {
                         const c = surfaceClass ? surfaceClass[r] : (e > 0 ? 1 : 0);
                         return c === 0 ? [0, 0, 0] : c === 2 ? [0.5, 0.5, 0.5] : [1, 1, 1];
                       }
                     : type === 'landmask' ? landMaskColor
                     : type === 'koppen'   ? (e, r) => koppenColor(koppen[r])
                     : type === 'biome'    ? (e, r) => biomeColor(koppen[r], e)
                     : (e) => elevationToColor(e);
            vertexColor = (s, it, ot, br, col) => {
                const c = fn(r_elevation[br], br);
                for (let k = 0; k < 3; k++) for (let ch = 0; ch < channels; ch++) col[k][ch] = c[ch];
            };
        }

        const tR = performance.now();
        process.stderr.write(args.quiet ? '' : `  rasterizing ${type} `);
        const img = rasterize(mesh, r_xyz, t_xyz, width, height, vertexColor, channels, background,
            spec.gray16, args.quiet ? null : () => process.stderr.write('.'));
        if (!args.quiet) process.stderr.write(` ${((performance.now() - tR) / 1000).toFixed(1)}s\n`);

        const file = path.join(outDir, `orogen-${spec.file}-${name}.png`);
        const tW = performance.now();
        let bytes;
        if (spec.gray16) {
            bytes = writePNG(file, width, height, 0, 16, (y, buf, off) => {
                const row = y * width;
                for (let x = 0; x < width; x++) {
                    const v = img[row + x];
                    buf[off + x * 2] = v >>> 8;
                    buf[off + x * 2 + 1] = v & 0xff;
                }
            });
        } else {
            bytes = writePNG(file, width, height, 6, 8, (y, buf, off) => {
                const row = y * width * 3;
                for (let x = 0; x < width; x++) {
                    const s = row + x * 3, d = off + x * 4;
                    buf[d] = img[s]; buf[d + 1] = img[s + 1]; buf[d + 2] = img[s + 2];
                    buf[d + 3] = 255;
                }
            });
        }
        log(`  wrote ${file} — ${humanBytes(bytes)} (encode ${((performance.now() - tW) / 1000).toFixed(1)}s)`);
        written.push(file);
    }

    // A planet code packs slider indices only — it cannot encode radius or
    // gravity, so on a non-Earth planet the filename can never fully identify
    // the world. This sidecar closes that gap.
    const manifest = {
        format: 'world-orogen-map-export',
        planetCode: effectiveCode,
        planetCodeGiven: args.code ?? null,
        planetCodeNote: 'A planet code encodes sliders only. Radius and gravity are NOT in it; '
                      + 'read them from planet below.',
        planet: planetSummary(ctx.planet),
        params: ctx.params,
        // Which drainage hypothesis this is a picture of. A basin list is not in
        // the planet code, so two maps of the same code can be two different
        // planets with nothing in the filename to say which. `finalElevation`
        // is the same hash the data export publishes, so a map can be tied to
        // the export it belongs with rather than assumed to match it.
        terrain: {
            finalElevation: hashTypedArray(ctx.r_elevation),
            basins: ctx.basins ? {
                selectionSource: ctx.basins.selectionSource,
                depressionsDetected: ctx.basins.catalogue.length,
                preserved: ctx.basins.selected.length,
                preserveListGiven: ctx.basins.drainageHypothesis.preserveListGiven,
                carveListGiven: ctx.basins.drainageHypothesis.carveListGiven,
                carvedByRetainZero: ctx.basins.drainageHypothesis.carvedByRetainZero,
                // The entries themselves are in the data export; here just their
                // identity, so a mismatched list is detectable.
                listHash: (args.preserveBasinList || args.carveBasinList)
                    ? hashJson(ctx.basins.drainageHypothesis.entries) : null,
                listFile: args.basinListFile,
            } : { selectionSource: 'off' },
        },
        width, height,
        projection: 'equirectangular (plate carrée), cell-centre registered',
        types: args.types,
        files: written.map(f => path.basename(f)),
        elevationNote: 'Heightmap PNGs encode PHYSICAL height via elevToHeightKm with this '
                     + "planet's 1/g relief scaling. Full-range: -5..6 km over 0..65535. "
                     + 'Land-only: 0..6 km over 0..65535, ocean black.',
    };
    const mf = path.join(outDir, 'manifest.json');
    fs.writeFileSync(mf, JSON.stringify(manifest, null, 2));
    log(`  wrote ${mf}`);

    console.log = stdout;
    if (args.quiet) for (const f of written) console.log(f);
}

main().catch(err => {
    console.error(`\nerror: ${err.message}`);
    if (process.env.DEBUG) console.error(err.stack);
    process.exit(1);
});
