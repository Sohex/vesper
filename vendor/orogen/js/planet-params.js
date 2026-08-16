// Planetary parameters.
//
// Orogen used to hardcode Earth: 6371 km, no gravity, no rotation. That is a
// problem when the next stage is ExoPlaSim, whose entire purpose is simulating
// arbitrary planets — the pipeline would have been modelling the climate of a
// super-Earth whose topography was built for Earth's gravity.
//
// WHAT ACTUALLY AFFECTS TERRAIN
//
//   radiusKm  — sets every physical length scale. Cell area, mean edge length,
//               basin area thresholds, scarp gradients, the km-per-hop
//               conversion for the distance fields. Threaded everywhere.
//   gravityMS2 — sets maximum relief. Crustal strength limits how high a
//               mountain can stand before it spreads under its own weight, and
//               to first order that ceiling goes as 1/g. Doubling gravity
//               roughly halves the mountains.
//
// WHAT DOES NOT
//
//   rotationPeriodHours, obliquityDeg, and the orbital elements are carried as
//   metadata and exported for ExoPlaSim to consume. Orogen does not use them:
//   they drive climate, and climate is frozen here. Do not wire them into
//   terrain — a rotation-dependent terrain model is not something this codebase
//   is trying to be.
//
// Earth defaults are exact, so a run that does not name a planet is
// bit-identical to one from before this module existed.

export const EARTH = {
    name: 'Earth',
    radiusKm: 6371,
    gravityMS2: 9.80665,
    rotationPeriodHours: 23.9344696,
    obliquityDeg: 23.4392811,
    eccentricity: 0.0167086,
    solarConstantWM2: 1361,
};

/** Reference gravity for the relief scaling below. */
export const REFERENCE_GRAVITY_MS2 = EARTH.gravityMS2;

/**
 * Build a validated planet definition. Anything unspecified falls back to Earth.
 *
 * @param {object} spec partial planet parameters
 * @returns {object} full planet definition plus derived quantities
 */
export function makePlanet(spec = {}) {
    const p = { ...EARTH, ...spec };

    const positive = (key, label) => {
        const v = Number(p[key]);
        if (!Number.isFinite(v) || v <= 0) {
            throw new Error(`Planet ${label} must be a positive number, got ${p[key]}`);
        }
        p[key] = v;
    };
    positive('radiusKm', 'radius');
    positive('gravityMS2', 'gravity');
    positive('rotationPeriodHours', 'rotation period');

    const finite = (key, label, lo, hi) => {
        const v = Number(p[key]);
        if (!Number.isFinite(v) || v < lo || v > hi) {
            throw new Error(`Planet ${label} must be between ${lo} and ${hi}, got ${p[key]}`);
        }
        p[key] = v;
    };
    finite('obliquityDeg', 'obliquity', 0, 180);
    finite('eccentricity', 'eccentricity', 0, 0.99);

    p.surfaceAreaKm2 = 4 * Math.PI * p.radiusKm * p.radiusKm;
    p.radiusEarth = p.radiusKm / EARTH.radiusKm;
    p.gravityEarth = p.gravityMS2 / EARTH.gravityMS2;

    // Maximum-relief multiplier. Crustal strength caps the load a mountain root
    // can support; that ceiling scales as 1/g, so a high-gravity world has
    // subdued topography and a low-gravity one has Olympus Mons.
    //
    // Applied to PHYSICAL height at the model-to-kilometres conversion, not to
    // the model's own dimensionless elevation parameter. Scaling the parameter
    // and then pushing it through the nonlinear height curve damps the effect
    // badly — a 4% parameter change came out as a 1.3% height change — and it
    // would perturb erosion, whose constants are calibrated in model units.
    // Land only: ocean depth is set by the isostatic balance of the water
    // column against the crust, where g cancels. Exactly 1 at Earth gravity.
    p.reliefScale = REFERENCE_GRAVITY_MS2 / p.gravityMS2;

    return p;
}

/** True when this planet is Earth in every respect terrain depends on. */
export function isEarthLike(planet) {
    return planet.radiusKm === EARTH.radiusKm && planet.gravityMS2 === EARTH.gravityMS2;
}

/**
 * Surface gravity from mass and radius, for callers who prefer to specify a
 * planet that way. mass in Earth masses, radius in Earth radii.
 */
export function gravityFromMassRadius(massEarth, radiusEarth) {
    return EARTH.gravityMS2 * massEarth / (radiusEarth * radiusEarth);
}

/** Compact summary for the export manifest. */
export function planetSummary(planet) {
    return {
        name: planet.name,
        radiusKm: planet.radiusKm,
        radiusEarth: planet.radiusEarth,
        gravityMS2: planet.gravityMS2,
        gravityEarth: planet.gravityEarth,
        surfaceAreaKm2: planet.surfaceAreaKm2,
        rotationPeriodHours: planet.rotationPeriodHours,
        obliquityDeg: planet.obliquityDeg,
        eccentricity: planet.eccentricity,
        solarConstantWM2: planet.solarConstantWM2,
        reliefScale: planet.reliefScale,
        usedByOrogen: ['radiusKm', 'gravityMS2'],
        passedThroughForDownstream: [
            'rotationPeriodHours', 'obliquityDeg', 'eccentricity', 'solarConstantWM2',
        ],
        note: 'radiusKm sets every length scale in the terrain model; gravityMS2 scales maximum '
            + 'relief as 1/g. The rest is metadata for ExoPlaSim — Orogen does not consume it, '
            + 'because it drives climate and climate is handled downstream.',
    };
}
