// Where the ice that carves this terrain is placed.
//
// This is the input to glacial erosion, and it is the only place the generator
// decides which cells carry ice. `terrain-post.js` consumes the index and does
// the carving; nothing else builds one.
//
// TWO SOURCES, AND ONLY ONE OF THEM IS A MEASUREMENT
//
//   'mask'      an ice index supplied by the caller, one value per mesh region,
//               in the same region order as the export's `raw/` arrays. This is
//               the physical path: the index comes from a climatology, so
//               obliquity, orbit, stellar spectrum and rotation are all in it
//               by construction rather than by proxy.
//
//   'heuristic' the generator's own latitude-and-elevation ramp. Earth-
//               calibrated, blind to temperature, and kept because the browser
//               app has no climatology to hand. It is a PLACEHOLDER, and the
//               function that builds it says so at the call site.
//
// WHAT THE HEURISTIC GETS RIGHT, WHICH IS NOT NOTHING
//
// The altitude gate is expressed in the model's dimensionless elevation
// parameter, and the obvious reading is that this makes the snowline track
// maximum relief, which goes as 1/g, rather than a freezing altitude. Write out
// both sides before believing it:
//
//   relief    physical height is `reliefScale * f(elev)` with
//             `reliefScale = g_ref/g`, because a crustal root fails at
//             sigma/(rho g). So the ceiling goes as 1/g.
//   freezing  the height at which a surface temperature reaches freezing is
//             (T_s - T_freeze)/Gamma, and a dry adiabat is Gamma = g/cp. So the
//             freezing height ALSO goes as 1/g, at fixed surface temperature
//             and composition.
//
// The two cancel. A dimensionless gate is therefore the gravity-invariant form,
// and re-anchoring it to a fixed physical altitude would introduce an error
// rather than remove one. That cancellation is exact only for a dry adiabat at
// Earth's cp and Earth's surface temperature, and this world's measured
// environmental lapse rate is shallower than g/cp implies, so a residual
// remains -- but it is a fraction of the gate, where the temperature term the
// heuristic is missing entirely is the whole of it.
//
// So the heuristic's defect is TEMPERATURE, not gravity, and the fix is the
// mask rather than a gravity term. `reliefScale` is accepted here so the
// reasoning above has somewhere to live and so a caller cannot conclude from
// the signature that gravity was never considered.
//
// NO TIME AXIS. An index says WHERE ice sits, never for how long. The strength
// slider that scales it is a declared choice, not a measurement; see
// `docs/src/reference/no-time-axis.md` in the consuming project.

import {
    GLACIAL_LAT_DIVISOR, GLACIAL_ELEV_LOW, GLACIAL_ELEV_HIGH,
    GLACIAL_ELEV_FACTOR_SCALE, GLACIAL_ELEV_FACTOR_LAT_BASE,
    GLACIAL_ELEV_FACTOR_LAT_SCALE,
} from './terrain-config.js';

function smoothstep(x, edge0, edge1) {
    const t = Math.max(0, Math.min(1, (x - edge0) / (edge1 - edge0)));
    return t * t * (3 - 2 * t);
}

/**
 * Check a supplied ice mask against the mesh it is about to be indexed into.
 *
 * The mask is matched to the terrain BY REGION INDEX, which is only meaningful
 * when both were built from the same mesh. A mesh is a pure function of the
 * seed and the region count, so those two are the identity that has to agree;
 * a mask from another build would otherwise apply one planet's ice to another
 * planet's mountains and produce an ordinary-looking result.
 *
 * @param {object} mask       { values: Float32Array|Array, numRegions, seed }
 * @param {number} numRegions the mesh being generated
 * @param {number} seed       the seed being generated
 * @returns {Float32Array} the values, as a typed array
 */
export function validateIceMask(mask, numRegions, seed) {
    if (!mask || !mask.values) throw new Error('Ice mask has no values');
    const values = mask.values instanceof Float32Array
        ? mask.values : Float32Array.from(mask.values);
    if (values.length !== numRegions) {
        throw new Error(`Ice mask has ${values.length} regions, the mesh has ${numRegions}. `
            + 'A mask is matched by region index and is only valid for the mesh it was built on.');
    }
    if (mask.numRegions !== undefined && mask.numRegions !== numRegions) {
        throw new Error(`Ice mask declares ${mask.numRegions} regions, the mesh has ${numRegions}`);
    }
    if (mask.seed !== undefined && mask.seed !== null && Number(mask.seed) !== Number(seed)) {
        throw new Error(`Ice mask was built for seed ${mask.seed}, this run is seed ${seed}. `
            + 'The mesh differs, so region indices do not refer to the same places.');
    }
    for (let r = 0; r < values.length; r++) {
        const v = values[r];
        if (!Number.isFinite(v) || v < 0 || v > 1) {
            throw new Error(`Ice mask value at region ${r} is ${v}; an ice index is in [0, 1]`);
        }
    }
    return values;
}

/**
 * Build the per-region ice index that glacial erosion carves with.
 *
 * @param {object} mesh
 * @param {Float32Array} r_xyz        unit-sphere positions
 * @param {Float32Array} r_elevation  model elevation parameter
 * @param {Uint8Array} r_isOcean
 * @param {number} glacialStrength    the declared slider, scales the index
 * @param {object} [opts]
 * @param {object} [opts.iceMask]     { values, numRegions, seed } — the physical path
 * @param {number} [opts.seed]        the run's seed, for the mask identity check
 * @param {number} [opts.reliefScale] g_ref/g; see the header on why it does not enter
 * @returns {{ idx: Float32Array, source: string }}
 */
export function buildGlacIdx(mesh, r_xyz, r_elevation, r_isOcean, glacialStrength, opts = {}) {
    const N = mesh.numRegions;
    const idx = new Float32Array(N);

    if (opts.iceMask) {
        const values = validateIceMask(opts.iceMask, N, opts.seed);
        for (let r = 0; r < N; r++) {
            if (r_isOcean[r]) continue;
            idx[r] = values[r] * glacialStrength;
        }
        return { idx, source: 'mask' };
    }

    // Earth-calibrated placeholder. At strength 1 glaciation starts at about 50
    // degrees latitude and at 0.5 about 70, which is the erosion slider setting
    // the ice line and is exactly why the mask path exists.
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
    return { idx, source: 'heuristic' };
}
