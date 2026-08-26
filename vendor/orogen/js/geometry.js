// Sphere geometry helpers shared by the simulation, the basin analysis and
// the exporter. DOM-free and dependency-free.
//
// COORDINATE CONVENTION. Geography is derived from the region position vector
// as lat = asin(y), lon = atan2(x, z). This is NOT the frame the Fibonacci
// point generator built internally (that one is z-up); everything
// user-visible and everything exported uses this y-up convention.

// Earth radius, kept as the default for every function below so a caller that
// names no planet behaves exactly as before. Non-Earth runs pass radiusKm
// explicitly — see js/planet-params.js.
export const PLANET_RADIUS_KM = 6371;

/** Total surface area of Earth, km². Use planet.surfaceAreaKm2 for other worlds. */
export const PLANET_AREA_KM2 = 4 * Math.PI * PLANET_RADIUS_KM * PLANET_RADIUS_KM;

const RAD2DEG = 180 / Math.PI;

/** Per-region latitude/longitude in degrees. */
export function regionLatLon(r_xyz) {
    const n = r_xyz.length / 3;
    const lat = new Float32Array(n);
    const lon = new Float32Array(n);
    for (let r = 0; r < n; r++) {
        const x = r_xyz[3 * r], y = r_xyz[3 * r + 1], z = r_xyz[3 * r + 2];
        lat[r] = Math.asin(Math.max(-1, Math.min(1, y))) * RAD2DEG;
        lon[r] = Math.atan2(x, z) * RAD2DEG;
    }
    return { lat, lon };
}

/**
 * Area of each dual (Voronoi-ish) cell in km².
 *
 * The cell boundary is the ring of triangle centres around the region. Area is
 * the sum of spherical triangles (region centre, vertex i, vertex i+1) using
 * the vector form of the spherical excess, so it is exact for the polygon
 * described by those vertices — the cell corners are triangle centroids rather
 * than circumcentres, which is the same approximation the renderer draws.
 */
export function regionCellArea(mesh, r_xyz, t_xyz, radiusKm = PLANET_RADIUS_KM) {
    const { numRegions, adjOffset } = mesh;
    const adjTri = mesh._adjTriList;
    const area = new Float32Array(numRegions);
    const R2 = radiusKm * radiusKm;

    const norm = (out, arr, i) => {
        let a = arr[3 * i], b = arr[3 * i + 1], c = arr[3 * i + 2];
        const l = Math.sqrt(a * a + b * b + c * c) || 1;
        out[0] = a / l; out[1] = b / l; out[2] = c / l;
        return out;
    };
    const p = [0, 0, 0], q = [0, 0, 0], s = [0, 0, 0];

    for (let r = 0; r < numRegions; r++) {
        norm(p, r_xyz, r);
        const start = adjOffset[r], end = adjOffset[r + 1];
        let total = 0;
        for (let i = start; i < end; i++) {
            const t0 = adjTri[i];
            const t1 = adjTri[i + 1 === end ? start : i + 1];
            norm(q, t_xyz, t0);
            norm(s, t_xyz, t1);
            // Spherical excess: E = 2·atan2(|p·(q×s)|, 1 + p·q + q·s + s·p)
            const cx = q[1] * s[2] - q[2] * s[1];
            const cy = q[2] * s[0] - q[0] * s[2];
            const cz = q[0] * s[1] - q[1] * s[0];
            const triple = Math.abs(p[0] * cx + p[1] * cy + p[2] * cz);
            const denom = 1 + (p[0] * q[0] + p[1] * q[1] + p[2] * q[2])
                            + (q[0] * s[0] + q[1] * s[1] + q[2] * s[2])
                            + (s[0] * p[0] + s[1] * p[1] + s[2] * p[2]);
            total += 2 * Math.atan2(triple, denom);
        }
        area[r] = total * R2;
    }
    return area;
}

/**
 * Mean spacing between adjacent region centres, km. The canonical resolution
 * scale used throughout the pipeline for converting cell-hop counts into
 * physical distance.
 */
export function avgEdgeKm(numRegions, radiusKm = PLANET_RADIUS_KM) {
    return (Math.PI * radiusKm) / Math.sqrt(numRegions);
}

/** Mean area of a single dual cell, km². */
export function avgCellAreaKm2(numRegions, radiusKm = PLANET_RADIUS_KM) {
    return (4 * Math.PI * radiusKm * radiusKm) / numRegions;
}

/**
 * Uniform per-cell area fallback, for callers that have no triangle centres.
 * Every cell gets the mesh mean — correct in aggregate, wrong per-cell.
 */
export function uniformCellArea(numRegions, radiusKm = PLANET_RADIUS_KM) {
    const a = avgCellAreaKm2(numRegions, radiusKm);
    const out = new Float32Array(numRegions);
    out.fill(a);
    return out;
}

// ─────────────────────────────────────────────────────────────────────────
//  Gaussian grids
// ─────────────────────────────────────────────────────────────────────────

/**
 * Gauss–Legendre latitudes, north to south, in degrees.
 *
 * Spectral models — ExoPlaSim among them — do not use equally spaced
 * latitudes. Their grid rows sit at the roots of the Legendre polynomial P_n,
 * which cluster differently from equal-angle spacing by a degree or more at mid
 * latitudes. Resampling equirectangular output onto that grid afterwards is an
 * extra lossy step; emitting it directly from the mesh avoids it.
 *
 * Roots are found by Newton iteration from the standard Chebyshev-like initial
 * guess, which converges in a handful of steps for any practical n.
 */
export function gaussianLatitudes(n) {
    const lat = new Float64Array(n);
    const weights = new Float64Array(n);
    const half = (n + 1) >> 1;

    for (let i = 0; i < half; i++) {
        // Initial guess for the i-th root of P_n, counting from +1.
        let z = Math.cos(Math.PI * (i + 0.75) / (n + 0.5));
        let dp = 0;
        for (let iter = 0; iter < 100; iter++) {
            // Legendre recurrence: p1 = P_n(z), p2 = P_{n-1}(z)
            let p1 = 1, p2 = 0;
            for (let j = 0; j < n; j++) {
                const p3 = p2;
                p2 = p1;
                p1 = ((2 * j + 1) * z * p2 - j * p3) / (j + 1);
            }
            dp = n * (z * p1 - p2) / (z * z - 1);
            const dz = p1 / dp;
            z -= dz;
            if (Math.abs(dz) < 1e-15) break;
        }
        const w = 2 / ((1 - z * z) * dp * dp);
        lat[i] = Math.asin(z) * RAD2DEG;              // northern half
        lat[n - 1 - i] = -lat[i];                     // symmetric southern half
        weights[i] = w;
        weights[n - 1 - i] = w;
    }
    return { lat, weights };
}

/**
 * Standard spectral truncations and their transform grids, as used by
 * PlaSim / ExoPlaSim. nlon is 3T+1 rounded up to an FFT-friendly size; nlat is
 * nlon/2.
 */
export const SPECTRAL_GRIDS = {
    T21:  { width: 64,  height: 32 },
    T31:  { width: 96,  height: 48 },
    T42:  { width: 128, height: 64 },
    T63:  { width: 192, height: 96 },
    T85:  { width: 256, height: 128 },
    T106: { width: 320, height: 160 },
    T127: { width: 384, height: 192 },
    T170: { width: 512, height: 256 },
};

/** Resolve "T42" (or "t42") to its grid, or null if it is not a truncation. */
export function spectralGrid(name) {
    if (typeof name !== 'string') return null;
    const key = name.trim().toUpperCase();
    return SPECTRAL_GRIDS[key] ? { ...SPECTRAL_GRIDS[key], truncation: key } : null;
}

/**
 * ASSIGNMENT boundaries for a latitude axis, as sine-of-latitude edges running
 * north to south. This is the nearest-row partition: a point belongs to the row
 * whose centre is closest to it in sin(lat). Binning by sin(lat) is what makes
 * assignment to a Gaussian row correct, because the rows are unequal in angle
 * but the edges are still monotonic.
 *
 * These are NOT the cell boundaries and their band widths are NOT cell areas.
 * On a Gaussian grid the cell is the quadrature interval, which is wider at the
 * pole than the nearest-row interval by 22 per cent at every truncation; on an
 * equally spaced grid the cell boundary is midway in LATITUDE, not in its sine.
 * `cellSinEdges` is the partition an area comes from. Two partitions, two
 * questions: which row does this point fall in, and how much sphere does this
 * row own.
 */
export function latitudeEdges(latDeg) {
    const n = latDeg.length;
    const edges = new Float64Array(n + 1);
    edges[0] = 1;                                     // north pole, sin(90) = 1
    for (let i = 1; i < n; i++) {
        edges[i] = (Math.sin(latDeg[i - 1] / RAD2DEG) + Math.sin(latDeg[i] / RAD2DEG)) / 2;
    }
    edges[n] = -1;                                    // south pole
    return edges;
}

/** Row index for a given latitude, given descending sine edges. Binary search. */
export function rowForLatitude(sinLat, edges) {
    let lo = 0, hi = edges.length - 1;                // edges descend from +1 to -1
    while (hi - lo > 1) {
        const mid = (lo + hi) >> 1;
        if (sinLat <= edges[mid]) lo = mid; else hi = mid;
    }
    return lo;
}

/**
 * CELL boundaries for a latitude axis, as sine-of-latitude edges running north
 * to south. This is the partition the rows OWN, and it is what an area and an
 * area weight come from.
 *
 * On a Gaussian grid the rows sit at quadrature abscissae rather than at cell
 * centres, so there is no geometric midpoint to appeal to and the partition is
 * chosen rather than derived. The one that matters is the quadrature's own: the
 * Gauss-Legendre weights ARE the sine extents of the intervals a spectral model
 * integrates over, they sum to 2, and laying them end to end from the north pole
 * gives the only partition on which a global mean of a model field equals the
 * model's own global mean.
 *
 * On an equally spaced grid the rows ARE cell centres, the boundary is midway in
 * latitude, and its sine is what the band area needs.
 */
export function cellSinEdges(grid) {
    const { height, type, lat, weights } = grid;
    const edges = new Float64Array(height + 1);
    edges[0] = 1;                                     // north pole, sin(90) = 1
    if (type === 'gaussian' && weights) {
        for (let j = 0; j < height; j++) edges[j + 1] = edges[j] - weights[j];
    } else {
        for (let j = 1; j < height; j++) {
            edges[j] = Math.sin((lat[j - 1] + lat[j]) / 2 / RAD2DEG);
        }
    }
    edges[height] = -1;                               // south pole
    return edges;
}

/**
 * Area of each cell of a lat/lon grid, km², row-major from the north pole.
 *
 * A cell spans a constant longitude width and a latitude band whose sine edges
 * `cellSinEdges` gives, and the area of such a band is R²·Δlon·Δ(sin lat) with no
 * approximation. Summed over the grid it is 4πR².
 *
 * This is the correct weight for area-averaging a gridded field, and on a
 * Gaussian grid it is `grid/gauss_weights.bin` times R²·Δlon, so the two agree
 * by construction rather than by luck.
 *
 * Three area quantities exist and they are three different things. The
 * per-region `cell_area` field is a MESH property: resampled onto a grid it
 * becomes the mean area of the mesh regions inside each cell, which is roughly
 * constant everywhere and carries no information about the cell, so it must NOT
 * be used to weight a gridded field. The nearest-row band widths in
 * `latitudeEdges` are an ASSIGNMENT partition and are not areas either.
 */
export function gridCellArea(grid, radiusKm = PLANET_RADIUS_KM) {
    const { width, height } = grid;
    const sinEdges = cellSinEdges(grid);
    const out = new Float32Array(width * height);
    const dLon = 2 * Math.PI / width;
    const R2 = radiusKm * radiusKm;
    for (let j = 0; j < height; j++) {
        const band = R2 * dLon * (sinEdges[j] - sinEdges[j + 1]);
        for (let i = 0; i < width; i++) out[j * width + i] = band;
    }
    return out;
}
