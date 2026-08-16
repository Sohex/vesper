// Lithology — surface rock type, as a two-layer model.
//
// WHAT THIS IS
//
// Rock type is a product of tectonic history, which is exactly what this repo
// simulates. Every input here already exists in the tectonics bundle: archetype
// weights, stress, subduction factor, boundary class, ridge/rift distances, and
// the hotspot / LIP / arc contributions. Nothing new is invented — this names
// what the model already decided, and hands the result to erosion so the
// landscape can respond to it.
//
// TWO LAYERS
//
//   basement — crystalline or oceanic crust, from deep tectonics
//   cover    — a sedimentary or volcanic veneer of some thickness on top
//
// Surface rock is the cover where it survives and the basement where it does
// not. Erosion strips cover; where thickness reaches zero, basement is exposed.
// That gives exhumed orogen cores and stripped cratons for free, without the
// cost of a real stratigraphic column — which this model could not support
// anyway, because Orogen is a snapshot with no time axis.
//
// WHAT THIS IS NOT
//
// It is not soil. Soil needs climate, vegetation and time (two of which arrive
// only after ExoPlaSim), so soil belongs to a later stage that combines this
// parent material with a real climate. It is also not geologic history: there
// are no ages, no stratigraphy and no unconformities, because there are no
// timesteps to hang them on. Seafloor age is the one exception and is a genuine
// ridge-distance proxy; continental age here is ordinal at best.

import {
    LITHO_CRATON_T, LITHO_FOLDBELT_T, LITHO_BASIN_T, LITHO_PLATEAU_T,
    LITHO_SUBDUCT_MELANGE_T, LITHO_ARC_T, LITHO_LIP_T, LITHO_HOTSPOT_T,
    LITHO_SHELF_DIST_CELLS, LITHO_CARBONATE_LAT_DEG,
    LITHO_COVER_BASIN_KM, LITHO_COVER_CRATON_KM, LITHO_COVER_PLATEAU_KM,
    LITHO_COVER_LIP_KM, LITHO_COVER_SHELF_KM, LITHO_COVER_PELAGIC_MAX_KM,
    LITHO_COVER_ARC_KM, LITHO_COVER_RIFT_KM, LITHO_COVER_EVAPORITE_KM,
    LITHO_SALT_CRUST_DEPTH_FRAC,
    LITHO_ERODIBILITY_STRENGTH,
    SCARP_MIN_GRADIENT, SCARP_FULL_GRADIENT, SCARP_EDGE_KM, SCARP_CONTRAST_SCALE,
} from './terrain-config.js';
import { avgEdgeKm, PLANET_RADIUS_KM } from './geometry.js';

/**
 * Rock classes.
 *
 * `erodibility` is a RELATIVE stream-power multiplier, not an absolute rate.
 * The field built from it is renormalised so its land mean is 1, so lithology
 * redistributes erosion rather than scaling it — the Hydraulic Erosion slider
 * keeps the meaning it always had, and the only new effect is differential
 * relief. Values are ordered by rock strength: unconsolidated sediment and
 * shale strip fast, quartzite and fresh granite hold up ridges.
 *
 * `densityGCm3` is carried for downstream use (isostasy, crustal budgets); the
 * pipeline itself does not consume it.
 *
 * `albedo` is bare-rock shortwave reflectance, a surface boundary condition for
 * ExoPlaSim. Vegetation and snow override it wherever they exist, so it matters
 * most on deserts, fresh lava and — strikingly — salt crust, which is among the
 * brightest natural land surfaces there are. Orogen does not apply it to
 * anything; it is exported for the climate stage to use.
 *
 * THAT MAKES SOME OF THESE NUMBERS LOAD-BEARING DOWNSTREAM. `evaporite` at 0.50
 * used to be applied to every cell of every closed basin, under the name
 * "Evaporite / playa fill" — a clean-halite reflectance spread over a surface
 * that is mostly not halite. On a planet with 12.65% of its land in closed
 * basins that single constant was worth ~1.7 W/m² per 0.10 of error, about twice
 * the entire radiative effect of exporting its lakes as dry ground. A closed
 * basin is zoned, so the class is now zoned too: see saltCrustMask().
 */
export const ROCK_CLASSES = [
    // id, code, name, category, erodibility, density
    { id:  0, code: 'water',        name: 'Open water (no exposed rock)',        category: 'none',        erodibility: 1.00, densityGCm3: null, albedo: 0.06 },
    { id:  1, code: 'morb',         name: 'Mid-ocean ridge basalt',              category: 'igneous',     erodibility: 0.70, densityGCm3: 2.9, albedo: 0.1 },
    { id:  2, code: 'oib',          name: 'Ocean island basalt',                 category: 'igneous',     erodibility: 0.80, densityGCm3: 2.9, albedo: 0.1 },
    { id:  3, code: 'flood_basalt', name: 'Continental flood basalt',            category: 'igneous',     erodibility: 0.65, densityGCm3: 2.9, albedo: 0.1 },
    { id:  4, code: 'arc_basalt',   name: 'Island-arc basalt / basaltic andesite', category: 'igneous',   erodibility: 0.85, densityGCm3: 2.8, albedo: 0.13 },
    { id:  5, code: 'arc_andesite', name: 'Continental-arc andesite / dacite',   category: 'igneous',     erodibility: 0.90, densityGCm3: 2.7, albedo: 0.2 },
    { id:  6, code: 'rift_bimodal', name: 'Rift bimodal volcanics',              category: 'igneous',     erodibility: 0.95, densityGCm3: 2.7, albedo: 0.16 },
    { id:  7, code: 'granite',      name: 'Granite',                             category: 'igneous',     erodibility: 0.40, densityGCm3: 2.65, albedo: 0.3 },
    { id:  8, code: 'granodiorite', name: 'Arc-root granodiorite',               category: 'igneous',     erodibility: 0.45, densityGCm3: 2.7, albedo: 0.28 },
    { id:  9, code: 'gneiss',       name: 'Cratonic gneiss',                     category: 'metamorphic', erodibility: 0.35, densityGCm3: 2.75, albedo: 0.28 },
    { id: 10, code: 'schist',       name: 'Orogenic schist / phyllite',          category: 'metamorphic', erodibility: 1.10, densityGCm3: 2.8, albedo: 0.22 },
    { id: 11, code: 'quartzite',    name: 'Quartzite',                           category: 'metamorphic', erodibility: 0.25, densityGCm3: 2.65, albedo: 0.35 },
    { id: 12, code: 'melange',      name: 'Subduction mélange / blueschist',     category: 'metamorphic', erodibility: 1.60, densityGCm3: 2.8, albedo: 0.18 },
    { id: 13, code: 'shelf_clastic', name: 'Shelf sandstone / shale',            category: 'sedimentary', erodibility: 2.20, densityGCm3: 2.5, albedo: 0.3 },
    { id: 14, code: 'carbonate',    name: 'Carbonate platform',                  category: 'sedimentary', erodibility: 1.30, densityGCm3: 2.7, albedo: 0.35 },
    { id: 15, code: 'foreland_clastic', name: 'Foreland molasse / flysch',       category: 'sedimentary', erodibility: 2.60, densityGCm3: 2.45, albedo: 0.28 },
    { id: 16, code: 'continental_clastic', name: 'Intracratonic clastics',       category: 'sedimentary', erodibility: 2.40, densityGCm3: 2.45, albedo: 0.28 },
    { id: 17, code: 'pelagic',      name: 'Pelagic ooze / abyssal clay',         category: 'sedimentary', erodibility: 3.00, densityGCm3: 2.0, albedo: 0.25 },
    { id: 18, code: 'evaporite',    name: 'Evaporite salt crust',                category: 'sedimentary', erodibility: 3.50, densityGCm3: 2.2, albedo: 0.5 },
    { id: 19, code: 'playa_clastic', name: 'Playa mud / alluvial fan fill',      category: 'sedimentary', erodibility: 2.80, densityGCm3: 2.1, albedo: 0.3 },
];

const BY_CODE = Object.fromEntries(ROCK_CLASSES.map(c => [c.code, c.id]));
const R = BY_CODE;

/** id → class record, for consumers reading an exported rock field. */
export const ROCK_BY_ID = ROCK_CLASSES.slice().sort((a, b) => a.id - b.id);

// ─────────────────────────────────────────────────────────────────────────
//  Cover sequence
// ─────────────────────────────────────────────────────────────────────────

/**
 * THE STRATIGRAPHIC ORDER OF THE COVER — youngest first. First match wins.
 *
 * This model has no time axis, so this table is the only place the question
 * "which deposit is on top?" is answered. It used to be answered by the order
 * the branches happened to be typed in, and that produced the same bug three
 * times in one file: an orogen exhuming the closed basin inside it, a closed
 * basin on oceanic crust surfacing as pelagic ooze, and a closed basin on a
 * flood-basalt province surfacing as basalt. Each was a chronology decided by
 * text layout.
 *
 * The table IS the chain — classifyLithology() walks it rather than duplicating
 * it. A parallel list that merely documents an if/else chain is the failure the
 * planet-code RADICES drift already taught us: two things that must agree,
 * maintained by hand, silently diverging.
 *
 * Adding a cover class means adding an entry here, which means deciding where in
 * the sequence it sits and what it `produces`. There is a test that fails if a
 * cover class reaches the surface without being declared.
 *
 *   domain    'ocean' | 'continent' | 'both' — which crust this can apply to
 *   when(c)   predicate on the per-cell context
 *   apply(c)  sets c.cover, c.thickKm, and c.depositional
 *   produces  every rock code this entry can assign, for the coverage test
 *
 * BASIN FILL IS FIRST, and that placement is the one load-bearing decision in
 * the table: a closed basin is defined by the condition that keeps operating
 * after every episodic process has stopped — no outlet, so whatever arrives
 * stays. Volcanism resurfaces a landscape once; a closed basin fills
 * continuously afterwards. Everything else here is ordinary parent material.
 */
export const COVER_SEQUENCE = [
    {
        code: 'basin_fill', domain: 'both',
        produces: ['evaporite', 'playa_clastic'],
        when: (c) => !!c.endorheic,
        apply: (c) => {
            // Closed basins are sediment and salt traps — no outlet means
            // everything delivered stays, and evaporation concentrates it. But
            // the trap is ZONED: dissolved load precipitates in the sump that
            // repeatedly floods and dries, while the margins take the clastic
            // load as playa mud and alluvial fans. Treating the whole basin as
            // salt crust overstated its albedo badly (see the note on
            // ROCK_CLASSES). saltCrust is built by saltCrustMask().
            c.cover = c.saltCrust ? BY_CODE.evaporite : BY_CODE.playa_clastic;
            c.thickKm = LITHO_COVER_EVAPORITE_KM;
            // Exempt from the fold-belt thinning below: this is where an
            // orogen's debris GOES.
            c.depositional = true;
        },
    },
    {
        code: 'hotspot_ocean', domain: 'ocean', produces: ['oib'],
        when: (c) => c.hotV > LITHO_HOTSPOT_T,
        apply: (c) => { c.cover = BY_CODE.oib; c.thickKm = LITHO_COVER_LIP_KM * 0.5; },
    },
    {
        code: 'lip_ocean', domain: 'ocean', produces: ['flood_basalt'],
        when: (c) => c.lipV > LITHO_LIP_T,
        apply: (c) => { c.cover = BY_CODE.flood_basalt; c.thickKm = LITHO_COVER_LIP_KM; },
    },
    {
        code: 'shelf', domain: 'ocean', produces: ['carbonate', 'shelf_clastic'],
        // Continental shelf: terrigenous input swamps pelagic settling.
        when: (c) => c.distCoast !== null && c.distCoast < c.shelfCells,
        apply: (c) => { c.cover = shelfClass(c.r_xyz, c.r); c.thickKm = LITHO_COVER_SHELF_KM; },
    },
    {
        code: 'pelagic', domain: 'ocean', produces: ['pelagic'],
        // Pelagic drape thickens away from the ridge: seafloor age grows with
        // ridge distance, and sediment accumulates with age. This is the one
        // genuine age signal in the model — and note that it is an age, not a
        // date. It orders nothing in this table.
        when: () => true,
        apply: (c) => { c.cover = BY_CODE.pelagic; c.thickKm = LITHO_COVER_PELAGIC_MAX_KM * c.age; },
    },
    {
        code: 'lip_continental', domain: 'continent', produces: ['flood_basalt'],
        when: (c) => c.lipV > LITHO_LIP_T,
        apply: (c) => { c.cover = BY_CODE.flood_basalt; c.thickKm = LITHO_COVER_LIP_KM; },
    },
    {
        code: 'arc', domain: 'continent', produces: ['arc_basalt', 'arc_andesite'],
        when: (c) => c.arcV > LITHO_ARC_T,
        apply: (c) => {
            c.cover = c.bothOcean ? BY_CODE.arc_basalt : BY_CODE.arc_andesite;
            c.thickKm = LITHO_COVER_ARC_KM;
        },
    },
    {
        code: 'rift', domain: 'continent', produces: ['rift_bimodal'],
        when: (c) => c.riftDist !== null && Number.isFinite(c.riftDist) && c.riftDist < c.riftHalfWidth,
        apply: (c) => { c.cover = BY_CODE.rift_bimodal; c.thickKm = LITHO_COVER_RIFT_KM; },
    },
    {
        code: 'basin_clastic', domain: 'continent',
        produces: ['foreland_clastic', 'continental_clastic'],
        when: (c) => c.basin > LITHO_BASIN_T,
        apply: (c) => {
            // Foreland basins sit against an orogen and fill with its debris;
            // intracratonic basins fill with quieter, more mature clastics.
            c.cover = c.foreland ? BY_CODE.foreland_clastic : BY_CODE.continental_clastic;
            c.thickKm = LITHO_COVER_BASIN_KM * c.basin;
        },
    },
    {
        code: 'coastal', domain: 'continent', produces: ['carbonate', 'shelf_clastic'],
        when: (c) => c.elev > 0 && c.distCoastLand !== null && c.distCoastLand < c.shelfCells,
        apply: (c) => { c.cover = shelfClass(c.r_xyz, c.r); c.thickKm = LITHO_COVER_SHELF_KM * 0.5; },
    },
    {
        code: 'plateau', domain: 'continent', produces: ['continental_clastic'],
        when: (c) => c.plateau > LITHO_PLATEAU_T,
        apply: (c) => { c.cover = BY_CODE.continental_clastic; c.thickKm = LITHO_COVER_PLATEAU_KM; },
    },
    {
        code: 'craton', domain: 'continent', produces: ['continental_clastic'],
        when: (c) => c.craton > LITHO_CRATON_T,
        apply: (c) => { c.cover = BY_CODE.continental_clastic; c.thickKm = LITHO_COVER_CRATON_KM * c.craton; },
    },
    {
        code: 'bare_basement', domain: 'continent', produces: [],
        // Nothing was deposited here: the basement IS the surface. `produces` is
        // empty because this entry assigns no class of its own — whatever the
        // basement chain decided stands.
        when: () => true,
        apply: (c) => { c.cover = c.basement; c.thickKm = 0; },
    },
];

// ─────────────────────────────────────────────────────────────────────────
//  Classification
// ─────────────────────────────────────────────────────────────────────────

/**
 * Assign basement rock, cover rock and cover thickness to every region.
 *
 * NON-MUTATING with respect to elevation. Reads the tectonics bundle, the
 * elevation-decomposition debug layers, and (optionally) the preserved-basin
 * state so closed basins can accumulate evaporites rather than ordinary clastics.
 *
 * @returns {{basement: Uint8Array, cover: Uint8Array, coverThicknessKm: Float32Array,
 *            surface: Uint8Array}}
 */
export function classifyLithology(mesh, r_xyz, r_elevation, tectonics, debugLayers, opts = {}) {
    const n = mesh.numRegions;
    const basins = opts.basins || null;
    const edgeKm = avgEdgeKm(n, opts.radiusKm);

    const basement = new Uint8Array(n);
    const cover = new Uint8Array(n);
    const coverThicknessKm = new Float32Array(n);

    const t = tectonics;
    const dl = debugLayers || {};
    const lip = dl.lip || null;
    const hotspot = dl.hotspot || null;
    const margins = dl.margins || null;
    const backArc = dl.backArc || null;

    // Stress normalised to [0,1] so the metamorphic grade threshold means the
    // same thing regardless of how hard this particular planet is colliding.
    const maxStress = t.maxStress > 0 ? t.maxStress : 1;

    // Shelf reach in cells, from a physical distance — a fixed hop count would
    // mean a different width at every resolution.
    const shelfCells = Math.max(1, Math.round(LITHO_SHELF_DIST_CELLS * 100 / edgeKm));

    const RAD2DEG = 180 / Math.PI;
    const endorheic = basins ? basins.r_isEndorheic : null;
    // Which of those cells is salt rather than clastic fill. Absent (no basin
    // pass, or an older caller) the whole basin falls back to clastic fill,
    // which is the conservative half of the pair — never silently the bright one.
    const saltCrust = basins ? basins.r_isSaltCrust : null;

    // One reused context object for the COVER_SEQUENCE walk. Reused rather than
    // allocated per cell because this loop runs 2.5M times on a full planet;
    // every field is overwritten before each walk, so nothing carries over.
    const c = {
        r: 0, r_xyz, shelfCells, basement: 0, elev: 0,
        craton: 0, basin: 0, plateau: 0, lipV: 0, hotV: 0, arcV: 0, age: 1,
        endorheic: 0, saltCrust: 0, bothOcean: 0,
        distCoast: null, distCoastLand: null, riftDist: null,
        riftHalfWidth: Number(t.riftHalfWidth) || 0, foreland: false,
        cover: 0, thickKm: 0, depositional: false,
    };

    for (let r = 0; r < n; r++) {
        const elev = r_elevation[r];
        const oceanic = t.r_isOcean[r] === 1;
        const craton = t.r_t_craton ? t.r_t_craton[r] : 0;
        const fold = t.r_t_foldBelt ? t.r_t_foldBelt[r] : 0;
        const basin = t.r_t_basin ? t.r_t_basin[r] : 0;
        const plateau = t.r_t_plateau ? t.r_t_plateau[r] : 0;
        const stress = t.r_stress[r] / maxStress;
        const subduct = t.r_subductFactor ? t.r_subductFactor[r] : 0;
        const bType = t.r_boundaryType ? t.r_boundaryType[r] : 0;
        const lipV = lip ? lip[r] : 0;
        const hotV = hotspot ? hotspot[r] : 0;
        const arcV = Math.max(margins ? margins[r] : 0, backArc ? backArc[r] : 0);

        let bm;
        // Basement: what is left when everything above is stripped away.
        if (oceanic) {
            bm = R.morb;
        } else {
            if (subduct > LITHO_SUBDUCT_MELANGE_T && bType === 1) bm = R.melange;
            else if (arcV > LITHO_ARC_T && t.r_hasOcean && t.r_hasOcean[r]) bm = R.granodiorite;
            else if (fold > LITHO_FOLDBELT_T) bm = stress > 0.6 ? R.gneiss : R.schist;
            else if (craton > LITHO_CRATON_T) bm = R.gneiss;
            else bm = R.granite;

            // Quartzite ridges: highly resistant metasediment surviving in the
            // cores of old fold belts. Rare by design — it is what holds up a
            // linear ridge long after the surrounding schist has gone.
            if (fold > LITHO_FOLDBELT_T && craton > LITHO_CRATON_T && stress > 0.5) bm = R.quartzite;
        }

        // Cover: walk COVER_SEQUENCE, youngest first, first match wins. The
        // order of that table is the model's entire answer to "which deposit is
        // on top", so it lives there and is walked here rather than being
        // restated as a chain of branches that could drift from it.
        //
        // ridgeHalfWidth is a SCALAR band width for the whole mesh, not a
        // per-region array, and ridgeDist is Infinity for cells no ridge
        // reached — which is the oldest crust, hence age 1.
        const ridgeD = t.ridgeDist ? t.ridgeDist[r] : Infinity;
        const ridgeHalf = Math.max(1, Number(t.ridgeHalfWidth) || 1);
        c.r = r;
        c.basement = bm;
        c.elev = elev;
        c.craton = craton; c.basin = basin; c.plateau = plateau;
        c.lipV = lipV; c.hotV = hotV; c.arcV = arcV;
        c.age = oceanic && Number.isFinite(ridgeD) ? Math.min(1, ridgeD / (ridgeHalf * 8)) : 1;
        c.endorheic = endorheic ? endorheic[r] : 0;
        c.saltCrust = saltCrust ? saltCrust[r] : 0;
        c.bothOcean = t.r_bothOcean ? t.r_bothOcean[r] : 0;
        c.distCoast = t.dist_coast ? t.dist_coast[r] : null;
        c.distCoastLand = t.dist_coast_land ? t.dist_coast_land[r] : null;
        c.riftDist = t.riftDist ? t.riftDist[r] : null;
        c.foreland = (t.coastConvergent && t.coastConvergent[r])
            || (t.dist_mountain && t.dist_mountain[r] < shelfCells * 2);
        c.cover = bm;
        c.thickKm = 0;
        c.depositional = false;

        const domain = oceanic ? 'ocean' : 'continent';
        for (let i = 0; i < COVER_SEQUENCE.length; i++) {
            const step = COVER_SEQUENCE[i];
            if (step.domain !== 'both' && step.domain !== domain) continue;
            if (!step.when(c)) continue;
            step.apply(c);
            break;
        }
        let cv = c.cover, thick = c.thickKm;
        const depositional = c.depositional;

        if (!oceanic) {
            // Active orogens are being exhumed, not buried: cover thins toward
            // zero as fold-belt intensity rises.
            //
            // EXCEPT in a closed basin, which is where the orogen's debris
            // GOES. The no-outlet condition that defines an endorheic basin
            // guarantees deposition — Altiplano, Qaidam, Tarim are intermontane
            // basins holding kilometres of fill inside the most actively
            // deforming belts on Earth. Thinning them was an oversight, not a
            // model: the rule and the basin branch were written together and
            // nobody asked what the rule did to the branch. It cost 7.7M km2
            // of a 2.5M-region Vesper — 15% of preserved-basin area, 2.4% of
            // land — which kept `cover_rock` = playa/evaporite but carried
            // thickness 0, so the surface read as orogenic basement. 80% of
            // that area had fold >= 1, where the multiplier is exactly zero.
            // It is worth 0.0027 of mean land rock albedo, which ExoPlaSim
            // integrates directly.
            if (fold > LITHO_FOLDBELT_T && !depositional) thick *= Math.max(0, 1 - fold);
        }

        basement[r] = bm;
        cover[r] = cv;
        // Defensive: every input above is a normalised weight or a guarded
        // distance, but a non-finite thickness would propagate through
        // erodibility into elevation, so it never leaves this loop.
        coverThicknessKm[r] = Number.isFinite(thick) && thick > 0 ? thick : 0;
    }

    const surface = surfaceRock(basement, cover, coverThicknessKm, r_elevation, mesh);
    return { basement, cover, coverThicknessKm, surface };
}

/**
 * Carbonate platform vs terrigenous clastics.
 *
 * The one classification here with a climate flavour: carbonate platforms are a
 * warm-water phenomenon. Latitude is pure geometry so this stays inside
 * Orogen's remit, but it is a proxy for sea-surface temperature, and a
 * downstream stage holding real ExoPlaSim SST can reclassify these cells with
 * better information. Flagged in the manifest for exactly that reason.
 */
function shelfClass(r_xyz, r) {
    const y = Math.max(-1, Math.min(1, r_xyz[3 * r + 1]));
    const latDeg = Math.abs(Math.asin(y) * 180 / Math.PI);
    return latDeg < LITHO_CARBONATE_LAT_DEG ? BY_CODE.carbonate : BY_CODE.shelf_clastic;
}

/** Exposed rock: cover where it survives, basement where it has been stripped. */
export function surfaceRock(basement, cover, coverThicknessKm, r_elevation, mesh) {
    const n = basement.length;
    const out = new Uint8Array(n);
    for (let r = 0; r < n; r++) {
        out[r] = coverThicknessKm[r] > 1e-4 ? cover[r] : basement[r];
    }
    return out;
}

// ─────────────────────────────────────────────────────────────────────────
//  Erodibility
// ─────────────────────────────────────────────────────────────────────────

/**
 * Per-region stream-power multiplier from the exposed rock.
 *
 * Renormalised so the area-weighted mean over land is exactly 1. Without that,
 * turning lithology on would change how much the planet erodes overall and
 * silently recalibrate the Hydraulic Erosion slider; with it, the total erosion
 * budget is preserved and lithology only decides *where* that budget is spent.
 *
 * LITHO_ERODIBILITY_STRENGTH blends between a uniform field (0) and the full
 * contrast implied by the rock table (1).
 */
export function buildErodibility(surface, r_elevation, cellArea, strength = LITHO_ERODIBILITY_STRENGTH) {
    const n = surface.length;
    const k = new Float32Array(n);
    const table = new Float32Array(ROCK_CLASSES.length);
    for (const c of ROCK_CLASSES) table[c.id] = c.erodibility;

    for (let r = 0; r < n; r++) k[r] = table[surface[r]] || 1;

    // Area-weighted mean over the EROSION DOMAIN — elevation > 0, which is what
    // the stream-power solve actually erodes. Deliberately not the subaerial
    // surface_class mask used for published fractions: this normalisation has to
    // match the set of cells the solve touches, and it is computed before any
    // basin is known. Same word, different set, on purpose.
    let wsum = 0, ksum = 0;
    for (let r = 0; r < n; r++) {
        if (r_elevation[r] <= 0) continue;
        const w = cellArea ? cellArea[r] : 1;
        wsum += w;
        ksum += w * k[r];
    }
    const mean = wsum > 0 ? ksum / wsum : 1;
    const inv = mean > 0 ? 1 / mean : 1;

    for (let r = 0; r < n; r++) {
        const normalised = k[r] * inv;
        k[r] = 1 + (normalised - 1) * strength;
        if (k[r] < 0.02) k[r] = 0.02;         // never fully un-erodible
    }
    return { erodibility: k, meanBeforeNormalisation: mean };
}

/**
 * Live lithology state handed to erodeComposite.
 *
 * The erosion loop updates cover thickness as material is removed and swaps a
 * cell to its basement erodibility once the cover is gone. That progressive
 * exhumation is the whole point of the two-layer model: it is what turns a
 * buried craton into a stripped shield and an orogen into an exposed core.
 * The per-iteration update is O(N) arithmetic, negligible beside the O(N log N)
 * sort the loop already performs.
 */
export function buildLithoState(litho, r_elevation, cellArea, strength = LITHO_ERODIBILITY_STRENGTH) {
    // Default here as well as in buildErodibility: an undefined strength would
    // otherwise reach the basement arithmetic below and produce NaN, which
    // propagates straight into elevation through the stream-power solve.
    if (!Number.isFinite(strength)) strength = LITHO_ERODIBILITY_STRENGTH;
    const n = litho.surface.length;
    const surfaceK = buildErodibility(litho.surface, r_elevation, cellArea, strength);

    // Basement erodibility on the same normalised scale, so a cell that strips
    // to basement mid-run does not jump onto a different calibration.
    const table = new Float32Array(ROCK_CLASSES.length);
    for (const c of ROCK_CLASSES) table[c.id] = c.erodibility;
    const mean = surfaceK.meanBeforeNormalisation;
    const inv = mean > 0 ? 1 / mean : 1;
    const basementK = new Float32Array(n);
    const coverK = new Float32Array(n);
    for (let r = 0; r < n; r++) {
        const nb = (table[litho.basement[r]] || 1) * inv;
        basementK[r] = Math.max(0.02, 1 + (nb - 1) * strength);
        const nc = (table[litho.cover[r]] || 1) * inv;
        coverK[r] = Math.max(0.02, 1 + (nc - 1) * strength);
    }

    return {
        erodibility: surfaceK.erodibility,
        basementErodibility: basementK,
        coverErodibility: coverK,
        coverThicknessKm: Float32Array.from(litho.coverThicknessKm),
        basement: litho.basement,
        cover: litho.cover,
        meanBeforeNormalisation: mean,
    };
}

/**
 * Strip cover by however much each cell eroded this iteration, and expose
 * basement where it runs out. Called by erodeComposite between passes.
 */
export function updateExhumation(state, elevBefore, r_elevation) {
    const n = state.coverThicknessKm.length;
    let exposed = 0;
    for (let r = 0; r < n; r++) {
        const removed = elevBefore[r] - r_elevation[r];
        if (removed <= 0) continue;
        const before = state.coverThicknessKm[r];
        if (before <= 0) continue;
        const after = before - removed;
        state.coverThicknessKm[r] = after > 0 ? after : 0;
        if (after <= 0) {
            state.erodibility[r] = state.basementErodibility[r];
            exposed++;
        }
    }
    return exposed;
}

/** Final surface rock after erosion has stripped what it stripped. */
export function finalSurfaceRock(state) {
    const n = state.coverThicknessKm.length;
    const out = new Uint8Array(n);
    for (let r = 0; r < n; r++) {
        out[r] = state.coverThicknessKm[r] > 1e-4 ? state.cover[r] : state.basement[r];
    }
    return out;
}

/**
 * Which cells of a preserved closed basin are salt crust rather than clastic fill.
 *
 * A closed basin is not a uniform surface. Runoff carries both a dissolved and a
 * clastic load; the clastic load drops at the basin margin as alluvial fans and
 * playa mud, while the dissolved load travels to the lowest ground and
 * precipitates there when the water evaporates. The salt crust is therefore the
 * SUMP — the part that repeatedly floods and dries — not the basin.
 *
 * Zoned by depth below the spill point rather than by area or cell count, so the
 * rule means the same thing at 2K regions and at 2.5M: a cell is salt crust when
 * it lies within LITHO_SALT_CRUST_DEPTH_FRAC of the basin's own relief above its
 * sink. A basin with no measurable relief is a pan that floods entirely, so all
 * of it is crust.
 *
 * This is deliberately a geometric proxy for flooding frequency and not a water
 * balance — how often a basin actually floods is climate, and climate is
 * downstream. What Orogen can say is where the low ground is.
 *
 * @param selected preserved basins, each with sink/spill elevations in model units
 * @param basinIndex per-region owning basin, -1 outside
 * @param r_elevation surface the basins were catalogued on
 */
export function saltCrustMask(selected, basinIndex, r_elevation, opts = {}) {
    const frac = Number.isFinite(opts.depthFraction) ? opts.depthFraction : LITHO_SALT_CRUST_DEPTH_FRAC;
    const n = r_elevation.length;
    const out = new Uint8Array(n);
    if (!selected || !basinIndex) return out;

    // Ceiling per basin, precomputed so the per-cell loop stays a comparison.
    const ceiling = new Float64Array(selected.length);
    for (let i = 0; i < selected.length; i++) {
        const b = selected[i];
        const sink = Number(b.sinkElevation);
        const depth = Number(b.depth);
        if (!Number.isFinite(sink)) { ceiling[i] = Infinity; continue; }
        ceiling[i] = Number.isFinite(depth) && depth > 0 ? sink + frac * depth : Infinity;
    }
    for (let r = 0; r < n; r++) {
        const i = basinIndex[r];
        if (i < 0 || i >= ceiling.length) continue;
        if (r_elevation[r] <= ceiling[i]) out[r] = 1;
    }
    return out;
}

/** Per-region bare-rock albedo from the exposed surface. */
export function buildAlbedo(surface) {
    const table = new Float32Array(ROCK_CLASSES.length);
    for (const c of ROCK_CLASSES) table[c.id] = c.albedo ?? 0.25;
    const out = new Float32Array(surface.length);
    for (let r = 0; r < surface.length; r++) out[r] = table[surface[r]] ?? 0.25;
    return out;
}

/** Area-weighted composition of the exposed surface, for the manifest. */
/**
 * Area-weighted rock-class composition over an explicit domain.
 *
 * @param isLand Uint8Array marking the cells to count, or null for the whole
 *   planet. NEVER derive this from `elevation > 0` here: that is `land_mask`,
 *   which drops the dry sub-sea-level closed-basin floors — 1.9% of a planet
 *   with preserved endorheic drainage — and the bias is not uniform. Evaporite
 *   is exactly the class that accumulates on those floors, so an elevation-sign
 *   denominator under-reported it by 2.2 points (18.6% vs 20.8%), which is the
 *   lithology most characteristic of the drainage regime being modelled. Pass
 *   the subaerial mask from surface_class.
 */
export function rockComposition(surface, isLand, cellArea) {
    const totals = new Float64Array(ROCK_CLASSES.length);
    let total = 0;
    for (let r = 0; r < surface.length; r++) {
        if (isLand && !isLand[r]) continue;
        const w = cellArea ? cellArea[r] : 1;
        totals[surface[r]] += w;
        total += w;
    }
    return ROCK_CLASSES
        .map(c => ({ code: c.code, name: c.name, category: c.category, areaKm2: totals[c.id],
                     fraction: total > 0 ? totals[c.id] / total : 0 }))
        .filter(x => x.areaKm2 > 0)
        .sort((a, b) => b.areaKm2 - a.areaKm2);
}

// ─────────────────────────────────────────────────────────────────────────
//  Scarp potential
// ─────────────────────────────────────────────────────────────────────────

/**
 * Where an escarpment would form — a marker, not geometry.
 *
 * Cuestas, mesas, plateau margins and stepped canyons are the landform class a
 * strength contrast between two layers produces, and they are the most legible
 * geological signature a planet has. Orogen cannot draw one: cells are tens of
 * kilometres across and a scarp is a sub-kilometre feature, so the cliff lives
 * entirely below the grid. Rather than fake it with a slope threshold the mesh
 * cannot represent (see the talus note in terrain-post.js), this reports where
 * the ingredients are, for a higher-resolution downstream pass or a renderer to
 * act on.
 *
 * Three factors, multiplied, each in [0,1]:
 *
 *   contrast — how different the cover and basement are in erodibility. No
 *              contrast, no differential retreat, no cliff.
 *   edge     — proximity to where the cover/basement interface daylights, i.e.
 *              where cover survives on one side and has been stripped on the
 *              other. A buried interface produces nothing at the surface.
 *   relief   — local steepness as a true gradient, so the measure means the same
 *              at every region count. Flat ground hosts no escarpment however
 *              strong the contrast.
 *
 * A NOTE ON WHAT KIND OF SCARP. In this model cover is always the weaker layer
 * (sedimentary and volcanic cover runs 0.65–3.5 on the erodibility scale;
 * crystalline basement runs 0.35–0.45). So these are stripped-edge scarps —
 * plateau margins and shield/cover boundaries of the Great Escarpment type —
 * rather than classic resistant-caprock cuestas, which would need a resistant
 * unit over a weak one and therefore a third layer. Treat the field as "an
 * escarpment belongs here", not as a claim about which way it faces.
 */
export function computeScarpPotential(mesh, r_elevation, lithoState, neighborDist, opts = {}) {
    const { numRegions, adjOffset, adjList } = mesh;
    // Subaerial, not elevation > 0: a dry closed-basin floor below sea level is
    // land, and a plateau margin standing inside one is still an escarpment.
    // opts.isLand comes from surface_class; the elevation-sign fallback is only
    // for callers with no basin pass.
    const isLand = opts.isLand ? (r) => opts.isLand[r] === 1 : (r) => r_elevation[r] > 0;
    const radiusKm = opts.radiusKm ?? PLANET_RADIUS_KM;
    const edgeKm = avgEdgeKm(numRegions, radiusKm);
    const out = new Float32Array(numRegions);
    if (!lithoState) return out;

    const cover = lithoState.coverThicknessKm;
    const coverK = lithoState.coverErodibility;
    const baseK = lithoState.basementErodibility;

    // ── edge: distance in hops from a cover/basement exposure boundary ──
    const exposed = new Uint8Array(numRegions);          // 1 = basement at surface
    for (let r = 0; r < numRegions; r++) exposed[r] = cover[r] > 1e-4 ? 0 : 1;

    const hops = Math.max(1, Math.round(SCARP_EDGE_KM / edgeKm));
    const dist = new Int32Array(numRegions).fill(-1);
    const queue = [];
    for (let r = 0; r < numRegions; r++) {
        if (!isLand(r)) continue;                        // scarps are a land feature
        for (let i = adjOffset[r], iEnd = adjOffset[r + 1]; i < iEnd; i++) {
            if (exposed[adjList[i]] !== exposed[r]) { dist[r] = 0; queue.push(r); break; }
        }
    }
    for (let qi = 0; qi < queue.length; qi++) {
        const r = queue[qi];
        if (dist[r] >= hops) continue;
        for (let i = adjOffset[r], iEnd = adjOffset[r + 1]; i < iEnd; i++) {
            const nb = adjList[i];
            if (dist[nb] >= 0 || !isLand(nb)) continue;
            dist[nb] = dist[r] + 1;
            queue.push(nb);
        }
    }

    const smoothstep = (x, a, b) => {
        const t = Math.max(0, Math.min(1, (x - a) / (b - a)));
        return t * t * (3 - 2 * t);
    };

    for (let r = 0; r < numRegions; r++) {
        if (!isLand(r) || dist[r] < 0) continue;

        const edge = 1 - dist[r] / (hops + 1);
        if (edge <= 0) continue;

        const contrast = Math.min(1, Math.abs(coverK[r] - baseK[r]) / SCARP_CONTRAST_SCALE);
        if (contrast <= 0) continue;

        // Steepest true gradient to any lower land neighbour. neighborDist is a
        // chord on the unit sphere, so multiplying by the planet radius turns
        // the ratio into a real rise/run.
        let steepest = 0;
        for (let i = adjOffset[r], iEnd = adjOffset[r + 1]; i < iEnd; i++) {
            const nb = adjList[i];
            // Land neighbours only. Measuring against a bathymetric neighbour
            // would score the continental slope as an escarpment — that drop is
            // the shelf/slope break, a different landform on a different scale,
            // and it would light up every coastline.
            if (r_elevation[nb] <= 0) continue;
            if (r_elevation[nb] >= r_elevation[r]) continue;
            const d = (neighborDist[i] || 1e-9) * radiusKm;
            const g = (r_elevation[r] - r_elevation[nb]) / d;
            if (g > steepest) steepest = g;
        }
        const relief = smoothstep(steepest, SCARP_MIN_GRADIENT, SCARP_FULL_GRADIENT);

        out[r] = contrast * edge * relief;
    }
    return out;
}
