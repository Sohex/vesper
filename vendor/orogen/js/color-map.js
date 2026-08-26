// Elevation → RGB colour mapping.

// Convert raw mesh elevation (nonlinear, 0-~1 for land) to physical height
// in kilometres.  Hybrid S-curve: quartic start gives extensive flatlands,
// steepest rise around t≈0.75, derivative→0 at top so peaks compress.
// Ocean (elev < 0) is mapped with a linear scale (~5 km at -0.5).
// Vertical scale for LAND that lies below sea level — dry closed-basin floors.
//
// The ocean branch below multiplies by 10 because the model's negative range is
// calibrated for abyssal bathymetry: -0.5 means -5 km of seabed. A continental
// closed basin is produced by the same machinery and inherits that range, but it
// is not seabed, and running it through the bathymetric scale put the deepest
// floor at -5.6 km — thirteen times the Dead Sea. At 1 km per model unit the
// same floor reads -562 m, which is the right order for a real endorheic basin
// (Dead Sea -430 m, Turfan -154 m, Qattara -133 m).
//
// This is a calibration choice, not a derivation: the model had no defined
// mapping for below-sea-level land, because before endorheic basins were
// preserved no such terrain survived to the export.
export const SUBSEA_LAND_KM_PER_UNIT = 1.0;

// The land curve's scale: the height it assigns to one full unit of the model's
// elevation parameter. The whole land branch is RELIEF_KM_PER_UNIT * s(elev),
// with s the dimensionless shape described below.
export const RELIEF_KM_PER_UNIT = 6;

/**
 * Model elevation → physical height in km.
 *
 * The land branch is RELIEF_KM_PER_UNIT * s(elev), where the shape function
 *
 *     s(t) = t⁴(5 − 4t),   s'(t) = 20t³(1 − t)
 *
 * is a Hermite interpolant DEFINED ON [0, 1]. s(0) = s'(0) = s''(0) = 0 gives
 * the extensive flatlands; s(1) = 1 with s'(1) = 0 gives the compressed peaks.
 * Both are endpoint conditions, so [0, 1] is the shape's domain and says
 * nothing about how much relief the terrain may have.
 *
 * The model's elevation parameter is NOT bounded by 1. `applyFinalShaping`
 * normalises land by RANK between its own min and max, peak compression is a
 * power rather than a clamp, and detail noise, warping, erosion and ridge
 * sharpening all run afterwards. Both registered builds carry land above 1, out
 * past 1.4 at the higher region count; the measurement is in
 * notes/audits/relief-curve-domain.md.
 *
 * s is not usable there. s'(t) is NEGATIVE above 1, so the polynomial turns
 * over: it falls back through 6 km, reaches zero at t = 1.25 and goes below sea
 * level beyond. The old `Math.min(elev, 1)` was a guard against that inversion,
 * not a physical ceiling on relief, and it published every cell above 1 at one
 * identical height — a plateau that downstream slope, drainage and hypsometry
 * all read as real, and which gave a closed basin with both its sink and its
 * spill above the clamp a depth of exactly zero.
 *
 * ABOVE THE DOMAIN THE SHAPE IS REPLACED BY ITS ARGUMENT, s(t) := t, so the
 * branch is RELIEF_KM_PER_UNIT * elev. That introduces no constant the curve
 * did not already have: it is continuous at t = 1, where s(1) = 1 already; its
 * gradient is the curve's own mean gradient across its domain, 6 km per model
 * unit; and it is strictly increasing for every t, so no two distinct
 * elevations share a height again.
 *
 * What it costs, stated because it is real: the gradient at the join steps from
 * 0 on the left, which is exactly the peak-compression condition s'(1) = 0, to
 * 6 km per model unit on the right. That kink is the same shape as the one this
 * function already carries at sea level, where the land branch arrives with
 * slope 0 and the ocean branch leaves with 10 km per unit, and it replaces a
 * gradient of exactly zero across the whole of the terrain it affects.
 *
 * The physical ceiling on land relief is a separate thing and lives in
 * `scaledHeightKm` below, as the 1/g factor out of sigma/(rho*g). It is not
 * this domain edge, and the two must not be conflated.
 *
 * @param elev   the model's dimensionless elevation parameter
 * @param isLand whether this cell is LAND. Defaults to the elevation-sign test,
 *               which is right everywhere except a dry closed-basin floor below
 *               sea level — pass surface_class === 1 to get those right. Getting
 *               it wrong reads a -560 m salt pan as -5.6 km of ocean trench.
 */
export function elevToHeightKm(elev, isLand) {
    const land = isLand === undefined ? elev > 0 : !!isLand;
    if (elev <= 0) {
        // Below sea level: seabed on the bathymetric scale, dry basin floor on
        // the land scale.
        return land ? elev * SUBSEA_LAND_KM_PER_UNIT : elev * 10;
    }
    if (!land) return elev * 10;      // shouldn't happen; keep it total
    if (elev >= 1) return RELIEF_KM_PER_UNIT * elev;   // 1.0→6, 1.25→7.5, 1.5→9
    const t2 = elev * elev;
    return RELIEF_KM_PER_UNIT * t2 * t2 * (5 - 4 * elev);  // 0→0, 0.25→0.09, 0.5→1.13, 0.75→3.80, 1.0→6
}

/**
 * Physical height in km, INCLUDING this planet's 1/g relief scaling.
 *
 * This is the single definition of a rule that was previously written out at
 * three separate call sites, and it is the one that looks like a bug if you meet
 * it without the reasoning:
 *
 *     positive heights scale by reliefScale; negative ones do not.
 *
 * That asymmetry is deliberate and physical, not an oversight.
 *
 *   LAND above sea level is limited by how much load a crustal root can carry
 *   before it fails. That ceiling is sigma/(rho*g), so it goes as 1/g: a
 *   high-gravity world has subdued topography, a low-gravity one gets Olympus
 *   Mons. Scaling is correct here.
 *
 *   OCEAN DEPTH is an isostatic balance between the water-plus-oceanic-crust
 *   column and the continental column. Write out that mass balance and g
 *   multiplies every term, so it CANCELS. Ocean depth is set by density
 *   contrasts and crustal thicknesses, and a planet at 2 g has the same ocean
 *   depth as one at 1 g. Scaling it would be wrong.
 *
 * The same cancellation is why a `reliefScale` of exactly 1 must be bit-identical
 * to omitting it: see tools/test-basins.mjs.
 *
 * A caveat that follows from taking the argument seriously, and is NOT modelled
 * here: isostatically compensated LAND -- a high plateau floating on thick crust
 * -- is gravity-independent for exactly the same reason the ocean is. Only the
 * strength-supported part of land relief should scale. Separating the two needs a
 * crustal-thickness field, which this model does not carry, so the uniform 1/g on
 * land is an upper bound on the correction and over-suppresses plateaus.
 *
 * @param elev        model elevation parameter
 * @param reliefScale REFERENCE_GRAVITY_MS2 / gravityMS2, exactly 1 for Earth
 * @param isLand      from surface_class, NOT the elevation sign; a dry
 *                    closed-basin floor below sea level is land
 */
export function scaledHeightKm(elev, reliefScale, isLand) {
    const h = elevToHeightKm(elev, isLand);
    return h > 0 ? h * reliefScale : h;
}

// Biome base colors indexed by Köppen class ID (satellite-view palette).
// 0=Ocean delegated, 1-30 = land biomes.
const BIOME_COLORS = [
    null,                        //  0 Ocean — handled separately
    [0.05, 0.30, 0.05],         //  1 Af   Tropical rainforest — deep emerald
    [0.08, 0.33, 0.07],         //  2 Am   Tropical monsoon — dense green
    [0.42, 0.50, 0.18],         //  3 Aw   Tropical savanna — yellow-green
    [0.82, 0.72, 0.50],         //  4 BWh  Hot desert — sandy tan
    [0.60, 0.55, 0.48],         //  5 BWk  Cold desert — gray-brown
    [0.72, 0.62, 0.30],         //  6 BSh  Hot steppe — dry gold
    [0.55, 0.52, 0.32],         //  7 BSk  Cold steppe — muted olive-tan
    [0.18, 0.42, 0.12],         //  8 Cfa  Humid subtropical — mid green
    [0.12, 0.38, 0.10],         //  9 Cfb  Oceanic — rich green
    [0.10, 0.28, 0.10],         // 10 Cfc  Subpolar oceanic — dark muted green
    [0.45, 0.48, 0.22],         // 11 Csa  Hot-summer Mediterranean — khaki-green
    [0.40, 0.45, 0.20],         // 12 Csb  Warm-summer Mediterranean — chaparral
    [0.35, 0.40, 0.20],         // 13 Csc  Cold-summer Mediterranean — darker khaki
    [0.20, 0.44, 0.14],         // 14 Cwa  Humid subtropical monsoon — mid green
    [0.15, 0.40, 0.12],         // 15 Cwb  Subtropical highland — green
    [0.12, 0.32, 0.10],         // 16 Cwc  Cold subtropical highland — dark green
    [0.12, 0.36, 0.08],         // 17 Dfa  Hot-summer continental — forest green
    [0.10, 0.32, 0.08],         // 18 Dfb  Warm-summer continental — forest green
    [0.06, 0.22, 0.08],         // 19 Dfc  Subarctic — dark spruce green
    [0.05, 0.18, 0.07],         // 20 Dfd  Extremely cold subarctic — very dark
    [0.38, 0.38, 0.18],         // 21 Dsa  Hot-summer continental dry — olive-brown
    [0.35, 0.35, 0.17],         // 22 Dsb  Warm-summer continental dry — olive-brown
    [0.08, 0.22, 0.08],         // 23 Dsc  Subarctic dry summer — dark green
    [0.06, 0.18, 0.07],         // 24 Dsd  Extremely cold subarctic dry — very dark
    [0.14, 0.36, 0.10],         // 25 Dwa  Hot-summer continental monsoon — forest green
    [0.12, 0.32, 0.09],         // 26 Dwb  Warm-summer continental monsoon
    [0.07, 0.22, 0.08],         // 27 Dwc  Subarctic monsoon — dark spruce
    [0.05, 0.18, 0.07],         // 28 Dwd  Extremely cold subarctic monsoon
    [0.35, 0.32, 0.22],         // 29 ET   Tundra — earthy brown (sparse moss/lichen on rock)
    [0.78, 0.80, 0.84],         // 30 EF   Ice cap — blue-tinted white
];

// Rocky/alpine mountain color for high-elevation blending.
const ROCK_COLOR = [0.42, 0.38, 0.32];

// Altitude thresholds (km) by Köppen group:
//   [alpine line, snow line]
// Alpine line: vegetation gives way to rocky alpine terrain.
// Snow line: permanent snow begins.
function altitudeThresholds(classId) {
    if (classId <= 0)  return [0, 0];           // Ocean
    if (classId <= 3)  return [3.5, 5.5];       // Tropical (A)
    if (classId <= 7)  return [3.0, 5.0];       // Arid (B)
    if (classId <= 16) return [2.0, 3.5];       // Temperate (C)
    if (classId <= 18 || classId === 21 || classId === 22 ||
        classId === 25 || classId === 26) return [1.5, 3.0];  // Continental humid (D*a, D*b)
    if (classId <= 28) return [0.8, 2.0];       // Subarctic (D*c, D*d)
    if (classId === 29) return [0.4, 1.5];      // Tundra (ET) — rocky higher up, snow only at peaks
    return [0, 0.5];                             // Ice cap (EF)
}

// Satellite-view biome color: realistic land colors based on Köppen class
// and elevation, with ocean delegated to the standard ocean palette.
export function biomeColor(koppenId, elevation) {
    // Ocean
    if (koppenId === 0 || elevation <= 0) return elevationToColor(elevation);

    const base = BIOME_COLORS[koppenId] || [0.30, 0.50, 0.20];
    const hKm = elevToHeightKm(elevation);
    const [alpineLine, snowLine] = altitudeThresholds(koppenId);

    let r = base[0], g = base[1], b = base[2];

    // Low-elevation subtle darkening for depth (0-200m)
    if (hKm < 0.2) {
        const dark = 0.93 + 0.07 * (hKm / 0.2);
        r *= dark; g *= dark; b *= dark;
    }

    // Mid-elevation: gentle darkening to show terrain relief (200m to alpine line)
    if (alpineLine > 0 && hKm > 0.2 && hKm < alpineLine) {
        const t = (hKm - 0.2) / (alpineLine - 0.2);
        const darken = 1.0 - t * 0.15; // up to 15% darker at alpine line
        r *= darken; g *= darken; b *= darken;
    }

    // Alpine zone: blend toward rocky brown-gray above the tree/vegetation line
    if (alpineLine > 0 && hKm > alpineLine) {
        const rockZone = snowLine > alpineLine ? snowLine - alpineLine : 2.0;
        const rockT = Math.min(1, (hKm - alpineLine) / rockZone);
        const s = rockT * rockT; // ease-in for gradual transition
        r = r + (ROCK_COLOR[0] - r) * s;
        g = g + (ROCK_COLOR[1] - g) * s;
        b = b + (ROCK_COLOR[2] - b) * s;
    }

    // Snow zone: blend toward white above the snow line
    if (snowLine > 0 && hKm > snowLine) {
        const snowT = Math.min(1, (hKm - snowLine) / 2.5);
        const s = snowT * snowT; // ease-in for gradual snow buildup
        r = r + (0.92 - r) * s;
        g = g + (0.93 - g) * s;
        b = b + (0.96 - b) * s;
    }

    return [r, g, b];
}

export function elevationToColor(e) {
    if (e < -0.50) return [0.04, 0.06, 0.30];
    if (e < -0.10) { const t=(e+0.50)/0.40; return [0.04+t*0.07,0.06+t*0.14,0.30+t*0.18]; }
    if (e <  0.00) { const t=(e+0.10)/0.10; return [0.11+t*0.19,0.20+t*0.22,0.48+t*0.12]; }
    if (e <  0.002){ const t=e/0.002;         return [0.72+t*0.08,0.68-t*0.02,0.46-t*0.10]; }
    if (e <  0.25) { const t=(e-0.002)/0.248; return [0.20-t*0.06,0.54-t*0.12,0.12+t*0.08]; }
    if (e <  0.50) { const t=(e-0.25)/0.25;  return [0.14+t*0.30,0.42-t*0.14,0.20-t*0.06]; }
    if (e <  0.75) { const t=(e-0.50)/0.25;  return [0.44+t*0.16,0.28+t*0.12,0.14+t*0.18]; }
    { const t=Math.min(1,(e-0.75)/0.20);      return [0.60+t*0.35,0.40+t*0.50,0.32+t*0.60]; }
}
