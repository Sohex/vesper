#!/usr/bin/env node
/**
 * Headless planet generation + full data export.
 *
 * Runs the same pipeline the browser worker runs (js/pipeline.js) and writes
 * every field the simulation produced — including all the tectonic internals
 * that used to stay inside assignElevation — as flat binaries with a JSON
 * manifest, and optionally as NetCDF for direct ingestion downstream.
 *
 * --grid repeats. Generation dominates the cost of an export and the grid has
 * no say in it, so one invocation generates the planet once and writes every
 * requested grid from it.
 *
 * Requires: npm i delaunator
 *
 * Usage:
 *   node tools/export-planet.mjs --out out/world-01
 *   node tools/export-planet.mjs --seed 12345 --regions 250000 --grid 1024x512 --netcdf
 *   node tools/export-planet.mjs --grid 128x64 --grid-method mean --only elevation,plates,tectonics
 *   node tools/export-planet.mjs --grid T42 --out a --grid T21 --out b --no-raw
 *   node tools/export-planet.mjs --list-fields
 *
 * Run with --help for the full option list.
 */

import fs from 'node:fs';
import path from 'node:path';
import Delaunator from 'delaunator';

import { setDelaunator } from '../js/sphere-mesh.js';
import { buildExportBundle, FIELD_GROUPS, gridCoords, zipStore } from '../js/data-export.js';
import { writeNetCDF, netcdfDtype } from './lib/netcdf-write.mjs';
import { ROCK_CLASSES } from '../js/lithology.js';
import { spectralGrid } from '../js/geometry.js';
import { decodePlanetCode } from '../js/planet-code.js';
import { basinAreaFromSlider, BASIN_MIN_DEPTH_KM, BASIN_MIN_AREA_KM2,
         BASIN_MIN_CELLS } from '../js/terrain-config.js';
import { parseBasinList } from '../js/basins.js';
import { gravityFromMassRadius, makePlanet, planetSummary } from '../js/planet-params.js';

setDelaunator(Delaunator);

// Defaults mirror the app's slider defaults.
const DEFAULTS = {
    regions: 100000,
    plates: 24,
    jitter: 0.75,
    noise: 1.0,
    continents: 5,
    continentSizeVariety: 0,
    landCoverage: 0.3,
    terrainWarp: 0,
    smoothing: 0,
    hydraulic: 0.4,
    thermal: 0.3,
    ridge: 0.3,
    glacial: 0.3,
    temperatureOffset: 0,
    precipitationOffset: 0,
    // TAKEN FROM terrain-config.js, NOT RESTATED. These three were literals
    // here, and basinMinAreaKm2 drifted to 850 in the module while this file
    // still said 1000. That is not merely a stale default: elevation.js reads
    // BASIN_MIN_AREA_KM2 and BASIN_MIN_CELLS directly, with no override, while
    // basins.js takes `opts.minAreaKm2 ?? BASIN_MIN_AREA_KM2` -- and this file
    // always passes the opt. So one generation decided which depressions the
    // mesh RESOLVES at 850 and which ones entered the CATALOGUE at 1000.
    // Importing them is what makes those two the same number by construction.
    basinMinDepthKm: BASIN_MIN_DEPTH_KM,
    basinMinAreaKm2: BASIN_MIN_AREA_KM2,
    basinMinCells: BASIN_MIN_CELLS,
};

const HELP = `
Headless World Orogen generation + data export.

  --out DIR              output directory (default: out/planet-<seed>).
                         One per --grid; see the note there.
  --code STR             planet code; supplies seed and every slider. Any explicit
                         flag given alongside it wins. Codes predating the basin and
                         rock-contrast sliders decode with both OFF, which is what
                         those planets were — pass the flags to enable them.
  --seed N               generation seed (default: random)
  --regions N            mesh region count (default: ${DEFAULTS.regions})
  --plates N             number of tectonic plates (default: ${DEFAULTS.plates})
  --continents N         number of continents (default: ${DEFAULTS.continents})
  --continent-variety F  continent size variety 0..1 (default: ${DEFAULTS.continentSizeVariety})
  --land-coverage F      target land fraction (default: ${DEFAULTS.landCoverage})
  --jitter F             mesh point jitter 0..1 (default: ${DEFAULTS.jitter})
  --noise F              tectonic noise magnitude (default: ${DEFAULTS.noise})
  --warp F               terrain warp strength (default: ${DEFAULTS.terrainWarp})
  --smoothing F          pre-erosion smoothing (default: ${DEFAULTS.smoothing})
  --hydraulic F          hydraulic erosion (default: ${DEFAULTS.hydraulic})
  --thermal F            thermal erosion (default: ${DEFAULTS.thermal})
  --ridge F              ridge sharpening (default: ${DEFAULTS.ridge})
  --glacial F            glacial erosion (default: ${DEFAULTS.glacial})
  --temp-offset F        climate temperature offset (default: 0)
  --precip-offset F      climate precipitation offset (default: 0)

  --no-basins            disable endorheic basin preservation (vanilla drainage)
  --preserve-basin ID    always preserve this basin id, bypassing the size floors
                         (repeatable; get ids from --list-basins)
  --basin-min-depth KM   depth-below-spill floor (default: ${DEFAULTS.basinMinDepthKm})
  --basin-min-area KM2   basin floor area floor (default: ${DEFAULTS.basinMinAreaKm2})
  --basin-min-cells N    resolution floor in mesh cells (default: ${DEFAULTS.basinMinCells})
  --water-level ID=KM    lake surface elevation for a preserved basin (repeatable).
                         Orogen never infers these; supply them or leave basins dry.
  --preserve-basins FILE explicit drainage hypothesis: the preserved set becomes
                         exactly the basins listed, ignoring the size floors. One id
                         per line, optional retain fraction 0-1 after it (1 = keep the
                         rim, 0 = carve it open, between = proportional incision).
                         '#' comments; a JSON array of ids or {id,retain} also works.
  --carve-basins FILE    same format, subtractive: threshold selection minus these.
  --list-basins          generate, print the basin catalogue, and exit

  --ice-mask FILE        where ice sits, one float32 per mesh region in region order,
                         with a FILE.json sidecar. Consumed at generation, so glacial
                         erosion carves what a climatology says is frozen instead of
                         the Earth-calibrated latitude ramp. Matched BY REGION INDEX,
                         and refused unless the sidecar's mesh identity matches.

  --no-lithology         uniform erodibility (upstream erosion, no rock types)
  --lithology-strength F 0 = uniform erodibility, 1 = full rock contrast (default: 1)
  --list-rocks           print the rock-class table and exit

  --planet NAME          label for the planet (default: Earth)
  --radius KM            planetary radius (default: 6371). Sets every length scale.
  --gravity MS2          surface gravity (default: 9.80665). Scales max relief as 1/g.
  --mass-radius M,R      surface gravity from Earth masses and Earth radii instead
  --rotation HOURS       rotation period, passed through for ExoPlaSim
  --obliquity DEG        axial tilt, passed through for ExoPlaSim
  --eccentricity E       orbital eccentricity, passed through for ExoPlaSim
  --solar-constant WM2   insolation, passed through for ExoPlaSim

  --grid WxH             equirectangular grid size (default: 1024x512)
  --grid T42             spectral truncation instead: T21 T31 T42 T63 T85 T106 T127 T170.
                         Emits Gauss-Legendre latitudes directly, skipping the lossy
                         equirect intermediate a spectral model would otherwise need.
                         REPEATABLE, and the reason to repeat it is cost: generation
                         dominates an export, so one invocation generates the planet
                         ONCE and writes every grid from it. Each --grid opens a target
                         and needs its own --out.
                         The EXPORT options -- --out, --raw, --no-raw, --no-grid,
                         --no-subgrid, --grid-method, --only, --netcdf, --zip -- bind to
                         the target the --grid before them opened; given before the
                         first --grid they set the default every target is created
                         from. Everything else describes the PLANET and is global, so
                         the targets of one invocation share a terrain hash by
                         construction.
  --no-subgrid           skip the sub-grid orography statistics
  --grid-method M        'mean' (area-weighted) or 'nearest' (default: mean)
  --no-grid              skip the gridded output
  --no-raw               skip the per-region output
  --raw                  keep it, against a --no-raw default (per target)
  --no-climate           skip wind/ocean/precipitation/temperature/Koppen entirely
  --only G1,G2           limit to field groups: ${Object.keys(FIELD_GROUPS).join(', ')}
  --netcdf               additionally write planet.nc (gridded fields, CF-ish)
  --zip                  write a single planet.zip instead of a directory tree
  --list-fields          print the field catalogue and exit (no generation)
  --quiet                suppress progress output
  --help
`;

/**
 * Read an ice mask: a flat float32 little-endian array, one value per mesh
 * region in region order, plus a JSON sidecar at the same path with `.json`
 * appended.
 *
 * The sidecar is not optional. A bare array of numbers cannot say which mesh it
 * was measured on, and a mask indexed into the wrong mesh applies one planet's
 * ice to another planet's mountains and returns an ordinary-looking result. The
 * sidecar carries `numRegions` and `seed`, which together are the mesh's whole
 * identity, and glacial-ice.js refuses on a mismatch rather than warning.
 *
 * Region order, never longitude: the export and every climatology label their
 * columns differently, and matching across that boundary by coordinate has
 * silently matched zero cells more than once.
 */
function readIceMask(file) {
    if (!file) throw new Error('--ice-mask expects a path');
    const sidecarPath = `${file}.json`;
    if (!fs.existsSync(sidecarPath)) {
        throw new Error(`--ice-mask needs a sidecar at ${sidecarPath}: a mask with no `
            + 'mesh identity cannot be checked against the mesh it is indexed into');
    }
    const meta = JSON.parse(fs.readFileSync(sidecarPath, 'utf8'));
    const buf = fs.readFileSync(file);
    if (buf.byteLength % 4 !== 0) {
        throw new Error(`--ice-mask ${file} is ${buf.byteLength} bytes, not a whole number of float32`);
    }
    const values = new Float32Array(buf.buffer, buf.byteOffset, buf.byteLength / 4);
    if (meta.numRegions !== undefined && meta.numRegions !== values.length) {
        throw new Error(`--ice-mask sidecar declares ${meta.numRegions} regions, `
            + `${file} holds ${values.length}`);
    }
    return {
        values,
        numRegions: meta.numRegions ?? values.length,
        seed: meta.seed ?? null,
        provenance: { file, ...meta },
    };
}

function parseArgs(argv) {
    // ONE export target: a grid, where it is written, and how. Everything
    // outside this object describes the planet, and the planet is generated once
    // however many targets there are.
    const newTarget = () => ({
        out: null,
        gridWidth: 1024, gridHeight: 512, gridType: 'uniform', gridTruncation: null,
        gridMethod: 'mean', grid: true, raw: true, subgrid: true,
        only: null, netcdf: false, zip: false,
    });
    const a = {
        seed: null, ...DEFAULTS,
        climate: true, listFields: false, quiet: false,
        preserveBasins: true, preserveBasin: [], waterLevels: {}, listBasins: false,
        lithology: true, lithologyStrength: undefined, listRocks: false,
        preserveBasinList: null, carveBasinList: null,
        planet: {}, iceMask: null,
        code: null, _explicit: new Set(),
        targets: [], targetDefaults: newTarget(),
    };
    // An export option belongs to the target the last --grid opened; before the
    // first --grid it sets the default later targets are created from. That is
    // the whole of the ordering rule, and it is what lets one invocation say
    // "T42 keeps raw/, the rest do not" without repeating the generation flags.
    const tgt = () => (a.targets.length ? a.targets[a.targets.length - 1] : a.targetDefaults);
    const num = (v, flag) => {
        const n = Number(v);
        if (!Number.isFinite(n)) throw new Error(`${flag} expects a number, got "${v}"`);
        return n;
    };
    for (let i = 2; i < argv.length; i++) {
        const f = argv[i];
        a._explicit.add(f);
        switch (f) {
            case '--help': case '-h': console.log(HELP); process.exit(0); break;
            case '--out': tgt().out = argv[++i]; break;
            case '--code': a.code = argv[++i]; break;
            case '--seed': a.seed = num(argv[++i], f); break;
            case '--regions': case '--n': a.regions = num(argv[++i], f); break;
            case '--plates': a.plates = num(argv[++i], f); break;
            case '--continents': a.continents = num(argv[++i], f); break;
            case '--continent-variety': a.continentSizeVariety = num(argv[++i], f); break;
            case '--land-coverage': a.landCoverage = num(argv[++i], f); break;
            case '--jitter': a.jitter = num(argv[++i], f); break;
            case '--noise': a.noise = num(argv[++i], f); break;
            case '--warp': a.terrainWarp = num(argv[++i], f); break;
            case '--smoothing': a.smoothing = num(argv[++i], f); break;
            case '--hydraulic': a.hydraulic = num(argv[++i], f); break;
            case '--thermal': a.thermal = num(argv[++i], f); break;
            case '--ridge': a.ridge = num(argv[++i], f); break;
            case '--glacial': a.glacial = num(argv[++i], f); break;
            case '--temp-offset': a.temperatureOffset = num(argv[++i], f); break;
            case '--precip-offset': a.precipitationOffset = num(argv[++i], f); break;
            case '--grid': {
                const spec = argv[++i];
                const t = { ...a.targetDefaults };
                const sp = spectralGrid(spec);
                if (sp) {
                    t.gridWidth = sp.width; t.gridHeight = sp.height;
                    t.gridType = 'gaussian'; t.gridTruncation = sp.truncation;
                } else {
                    const m = /^(\d+)x(\d+)$/.exec(spec);
                    if (!m) throw new Error(`--grid expects WxH or a truncation like T42, got "${spec}"`);
                    t.gridWidth = +m[1]; t.gridHeight = +m[2];
                    t.gridType = 'uniform'; t.gridTruncation = null;
                }
                a.targets.push(t);
                break;
            }
            case '--no-subgrid': tgt().subgrid = false; break;
            case '--planet': a.planet.name = argv[++i]; break;
            case '--radius': a.planet.radiusKm = num(argv[++i], f); break;
            case '--gravity': a.planet.gravityMS2 = num(argv[++i], f); break;
            case '--mass-radius': {
                const parts = String(argv[++i]).split(',');
                if (parts.length !== 2) throw new Error('--mass-radius expects M,R in Earth units');
                a.planet.gravityMS2 = gravityFromMassRadius(Number(parts[0]), Number(parts[1]));
                break;
            }
            case '--rotation': a.planet.rotationPeriodHours = num(argv[++i], f); break;
            case '--obliquity': a.planet.obliquityDeg = num(argv[++i], f); break;
            case '--eccentricity': a.planet.eccentricity = num(argv[++i], f); break;
            case '--solar-constant': a.planet.solarConstantWM2 = num(argv[++i], f); break;
            case '--grid-method': tgt().gridMethod = argv[++i]; break;
            case '--no-grid': tgt().grid = false; break;
            case '--no-raw': tgt().raw = false; break;
            case '--raw': tgt().raw = true; break;
            case '--no-climate': a.climate = false; break;
            case '--only': tgt().only = argv[++i].split(',').map(s => s.trim()).filter(Boolean); break;
            case '--netcdf': tgt().netcdf = true; break;
            case '--zip': tgt().zip = true; break;
            case '--list-fields': a.listFields = true; break;
            case '--no-basins': a.preserveBasins = false; break;
            case '--list-basins': a.listBasins = true; break;
            case '--preserve-basins': a.preserveBasinList = parseBasinList(fs.readFileSync(argv[++i], 'utf8')); break;
            case '--carve-basins': a.carveBasinList = parseBasinList(fs.readFileSync(argv[++i], 'utf8')); break;
            case '--ice-mask': a.iceMask = readIceMask(argv[++i]); break;
            case '--no-lithology': a.lithology = false; break;
            case '--lithology-strength': a.lithologyStrength = num(argv[++i], f); break;
            case '--list-rocks': a.listRocks = true; break;
            case '--preserve-basin': a.preserveBasin.push(argv[++i]); break;
            case '--basin-min-depth': a.basinMinDepthKm = num(argv[++i], f); break;
            case '--basin-min-area': a.basinMinAreaKm2 = num(argv[++i], f); break;
            case '--basin-min-cells': a.basinMinCells = num(argv[++i], f); break;
            case '--water-level': {
                const spec = argv[++i];
                const eq = spec ? spec.lastIndexOf('=') : -1;
                if (eq <= 0) throw new Error(`--water-level expects ID=KM, got "${spec}"`);
                const id = spec.slice(0, eq);
                const km = Number(spec.slice(eq + 1));
                if (!Number.isFinite(km)) throw new Error(`--water-level: "${spec.slice(eq + 1)}" is not a number`);
                a.waterLevels[id] = km;
                break;
            }
            case '--quiet': a.quiet = true; break;
            default: throw new Error(`Unknown option: ${f} (try --help)`);
        }
    }
    // No --grid at all is one target on the default grid, which is what this
    // tool has always done.
    if (!a.targets.length) a.targets.push(a.targetDefaults);
    const byOut = new Map();
    for (const t of a.targets) {
        const where = a.targets.length > 1 ? ` for --grid ${gridLabel(t)}` : '';
        if (t.gridMethod !== 'mean' && t.gridMethod !== 'nearest') {
            throw new Error(`--grid-method must be 'mean' or 'nearest'`);
        }
        if (t.only) {
            const bad = t.only.filter(g => !FIELD_GROUPS[g]);
            if (bad.length) throw new Error(`Unknown field group(s): ${bad.join(', ')}. Valid: ${Object.keys(FIELD_GROUPS).join(', ')}`);
        }
        if (!t.raw && !t.grid) throw new Error(`Nothing to export${where}: --no-raw and --no-grid together`);
        if (a.targets.length === 1) continue;
        // One directory per grid. Letting two targets share a path would have
        // them overwrite each other field by field and leave a manifest
        // describing whichever ran last: a wrong export, not a failed one.
        if (!t.out) {
            throw new Error(`--grid ${gridLabel(t)} has no --out. With more than one --grid `
                + 'each target needs its own output directory.');
        }
        const key = path.resolve(t.out);
        if (byOut.has(key)) {
            throw new Error(`--grid ${gridLabel(t)} and --grid ${byOut.get(key)} both write to `
                + `${t.out}. One output directory per grid.`);
        }
        byOut.set(key, gridLabel(t));
    }
    return a;
}

/** How a target names itself in a message: its truncation, or WxH. */
function gridLabel(t) {
    return t.gridTruncation ?? `${t.gridWidth}x${t.gridHeight}`;
}

/** A NetCDF variable name has to be a valid identifier. */
function ncName(name) {
    const s = name.replace(/[^A-Za-z0-9_]/g, '_');
    return /^[A-Za-z_]/.test(s) ? s : `v_${s}`;
}

function writeNetCDFFile(outDir, manifest, target, files) {
    const { gridWidth: W, gridHeight: H } = target;
    const coords = gridCoords(W, H, target.gridType);
    const byPath = new Map(files.map(f => [f.name, f.bytes]));

    const variables = [
        { name: 'lat', dtype: 'float64', dimensions: ['lat'], data: coords.lat,
          attributes: { units: 'degrees_north', long_name: 'latitude', standard_name: 'latitude', axis: 'Y' } },
        { name: 'lon', dtype: 'float64', dimensions: ['lon'], data: coords.lon,
          attributes: { units: 'degrees_east', long_name: 'longitude', standard_name: 'longitude', axis: 'X' } },
    ];

    for (const entry of manifest.grid.fields) {
        if (entry.name === 'lat' || entry.name === 'lon') continue;
        const bytes = byPath.get(entry.path);
        if (!bytes) continue;
        const dtype = netcdfDtype(entry.dtype);
        const data = typedFromBytes(bytes, entry.dtype);
        variables.push({
            name: ncName(entry.name),
            dtype,
            dimensions: ['lat', 'lon'],
            data,
            attributes: {
                units: entry.units === '1' ? '1' : entry.units,
                long_name: entry.description,
                orogen_group: entry.group,
                orogen_resample: entry.method,
            },
        });
    }

    const nc = writeNetCDF({
        dimensions: { lat: H, lon: W },
        variables,
        attributes: {
            Conventions: 'CF-1.8',
            title: `World Orogen planet ${manifest.seed}`,
            source: 'World Orogen (worldbuilding fork), tools/export-planet.mjs',
            seed: String(manifest.seed),
            num_regions: String(manifest.numRegions),
            planet_radius_km: manifest.planetRadiusKm,
            comment: 'Equirectangular, cell-centre registered. Elevation in km, sea level at 0. '
                   + 'Categorical fields are nearest-region sampled; continuous fields are area-weighted means.',
        },
    });
    fs.writeFileSync(path.join(outDir, 'planet.nc'), nc);
    return nc.length;
}

function typedFromBytes(bytes, dtype) {
    const buf = bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength);
    switch (dtype) {
        case 'float64': return new Float64Array(buf);
        case 'int32': return new Int32Array(buf);
        case 'uint32': return new Uint32Array(buf);
        case 'int16': return new Int16Array(buf);
        case 'uint16': return new Uint16Array(buf);
        case 'int8': return new Int8Array(buf);
        case 'uint8': return new Uint8Array(buf);
        default: return new Float32Array(buf);
    }
}

function humanBytes(n) {
    const u = ['B', 'KB', 'MB', 'GB'];
    let i = 0;
    while (n >= 1024 && i < u.length - 1) { n /= 1024; i++; }
    return `${n.toFixed(i === 0 ? 0 : 1)} ${u[i]}`;
}

function printBasins(b, out = console.log) {
    const pad = (v, n) => String(v).padStart(n);
    const console = { log: out };   // shadow so --quiet cannot swallow the listing
    console.log(`\n${b.catalogue.length} depressions detected; ${b.selected.length} preserved.`);
    console.log(`Selection floors: depth >= ${b.resolution.minDepthKm} km, area >= ${b.resolution.minAreaKm2} km2, `
              + `>= ${b.resolution.minCells} cells (mean cell ${b.resolution.avgCellAreaKm2.toFixed(0)} km2).`);
    console.log('\nPreserved basins, largest first by volume:\n');
    console.log('  id                              depth_m    area_km2  catch_km2  retain  spills_into');
    for (let i = 0; i < b.selected.length; i++) {
        const s = b.selected[i], f = b.finalMeasurements[i], c = b.finalCatchments[i];
        console.log(`  ${s.id.padEnd(30)} ${pad((s.depthKm * 1000).toFixed(0), 7)} `
            + `${pad(Math.round(s.areaKm2), 11)} ${pad(Math.round(c.catchmentAreaKm2), 10)} `
            + `${pad((100 * (f.retainedFraction ?? 0)).toFixed(0) + '%', 7)}  ${s.spillsInto ?? '-'}`);
    }
    console.log('\nRetention is the fraction of the basin\'s natural spill depth surviving erosion.');
    console.log('Pass an id to --preserve-basin to force-keep one below the size floors,');
    console.log('or to --water-level ID=KM to fill it. Orogen does not decide lake levels.');
}

function listRocks() {
    console.log('\nRock classes. `erodibility` is a relative stream-power multiplier; the field');
    console.log('built from it is mean-normalised to 1 over land, so lithology decides where');
    console.log('erosion goes, not how much of it there is.\n');
    console.log('  id  code                  erod  dens  category      name');
    for (const c of ROCK_CLASSES) {
        console.log(`  ${String(c.id).padStart(2)}  ${c.code.padEnd(20)} ${c.erodibility.toFixed(2).padStart(5)} `
            + `${(c.densityGCm3 ?? '-').toString().padStart(5)}  ${c.category.padEnd(13)} ${c.name}`);
    }
    console.log('\nSurface rock is the cover where it survives and the basement where erosion');
    console.log('has stripped it. This is parent material, not soil — soil needs climate.');
}

function listFields() {
    for (const [group, names] of Object.entries(FIELD_GROUPS)) {
        console.log(`\n${group}`);
        for (const n of names) console.log(`  ${n}`);
    }
    console.log('\nFields present in a given run depend on the parameters — climate fields');
    console.log('need climate enabled, super-plate fields need --plates >= 8. Anything the');
    console.log('pipeline produces but this list omits is still exported, under group "other".');
}

async function main() {
    const args = parseArgs(process.argv);
    if (args.listFields) { listFields(); return; }
    if (args.listRocks) { listRocks(); return; }

    // A planet code fills in the sliders; anything given explicitly overrides it.
    if (args.code) {
        const d = decodePlanetCode(args.code);
        if (!d) throw new Error(`Could not decode planet code "${args.code}"`);
        const set = (flag, key, value) => { if (!args._explicit.has(flag)) args[key] = value; };
        set('--seed', 'seed', d.seed);
        set('--regions', 'regions', d.N);
        set('--plates', 'plates', d.P);
        set('--continents', 'continents', d.numContinents);
        set('--continent-variety', 'continentSizeVariety', d.continentSizeVariety);
        set('--land-coverage', 'landCoverage', d.landCoverage);
        set('--jitter', 'jitter', d.jitter);
        set('--noise', 'noise', d.roughness);
        set('--warp', 'terrainWarp', d.terrainWarp);
        set('--smoothing', 'smoothing', d.smoothing);
        set('--hydraulic', 'hydraulic', d.hydraulicErosion);
        set('--thermal', 'thermal', d.thermalErosion);
        set('--ridge', 'ridge', d.ridgeSharpening);
        set('--glacial', 'glacial', d.glacialErosion);
        set('--temp-offset', 'temperatureOffset', d.temperatureOffset);
        set('--precip-offset', 'precipitationOffset', d.precipitationOffset);
        args.toggledIndices = d.toggledIndices || [];
        // Any explicit basin or lithology flag means the caller is driving that
        // feature and the code must not override it. Listing only ONE of each
        // group is a trap: passing --basin-min-cells alongside a legacy code
        // (whose slider says "off") silently disabled basin preservation
        // entirely while looking like it had configured it.
        const BASIN_FLAGS = ['--no-basins', '--basin-min-area', '--basin-min-cells',
                             '--basin-min-depth', '--preserve-basin',
                             '--preserve-basins', '--carve-basins'];
        const LITHO_FLAGS = ['--no-lithology', '--lithology-strength'];
        const touched = (flags) => flags.some(f => args._explicit.has(f));

        if (!touched(BASIN_FLAGS)) {
            const basinArea = basinAreaFromSlider(d.basinSlider ?? 0);
            // THE CODE'S SLIDER MAY NOT SILENTLY OUTRANK THE FORK'S FLOOR.
            // basinAreaFromSlider returns a rung of BASIN_AREA_LADDER_KM2, and
            // BASIN_MIN_AREA_KM2 is a decision this fork took that need not be
            // on that ladder -- 850 is not. Taking the rung anyway splits the
            // generation in two, because elevation.js reads the constant
            // DIRECTLY for its resolvable test while the catalogue is selected
            // at whatever lands here. That ran twice before it was caught, at
            // 23 minutes an export, and the only evidence was minAreaKm2 in the
            // manifest. So say so and stop, naming the flag that settles it.
            if (basinArea > 0 && basinArea !== BASIN_MIN_AREA_KM2) {
                console.error(
                    `The planet code's basin slider gives a floor of ${basinArea} km2 and\n` +
                    `terrain-config.js declares BASIN_MIN_AREA_KM2 = ${BASIN_MIN_AREA_KM2}.\n` +
                    `js/elevation.js reads the constant directly, so taking the code's value\n` +
                    `here would select the catalogue at one floor and decide what the mesh\n` +
                    `resolves at another, in the same generation.\n\n` +
                    `Pass --basin-min-area explicitly to say which you mean. The fork's\n` +
                    `answer is --basin-min-area ${BASIN_MIN_AREA_KM2}; the code's own is\n` +
                    `--basin-min-area ${basinArea}.`);
                process.exit(2);
            }
            args.preserveBasins = basinArea > 0;
            if (basinArea > 0) args.basinMinAreaKm2 = basinArea;
        }
        if (!touched(LITHO_FLAGS)) {
            const strength = d.lithologyStrength ?? 0;
            args.lithology = strength > 0;
            args.lithologyStrength = strength;
        }
    }

    const seed = args.seed ?? Math.floor(Math.random() * 16777216);
    const defaultOut = path.join('out', `planet-${args.code || seed}`);

    // AN EXISTING EXPORT IS NEVER OVERWRITTEN, and the refusal is here rather
    // than at the write because generation is almost all of what an export
    // costs and a target that already holds a manifest is knowable before the
    // first plate moves. The basin-floor refusal above is here for the same
    // reason and cost two 23-minute exports to learn it.
    //
    // WHY REFUSING BY EXISTENCE IS THE RIGHT INSTRUMENT. An export directory is
    // what a world is RECONSTRUCTED from; the payload is large enough that no
    // consumer keeps it under version control, so a directory rewritten in
    // place has nothing to fall back on. Worse, the rewrite need not announce
    // itself: re-exporting the same code at the same seed onto a different grid
    // leaves a directory whose terrain hash still matches every reader's
    // record while its bytes have moved, and a reader that checks the hash
    // learns nothing. The write is the only moment the difference is visible.
    //
    // NO OVERRIDE FLAG, deliberately. A flag makes the unrecoverable write the
    // shorter path, and it buys nothing that removing the directory by hand
    // does not: a failed export is a directory whose caller deletes it, having
    // looked at what is in it. --list-basins is exempt because it generates and
    // prints without writing a target at all.
    if (!args.listBasins) {
        for (const target of args.targets) {
            const outDir = target.out ?? defaultOut;
            const held = ['manifest.json', 'planet.zip']
                .filter(name => fs.existsSync(path.join(outDir, name)));
            if (!held.length) continue;
            console.error(
                `${outDir} already holds an export (${held.join(', ')}).\n` +
                `Exporting into it would replace a build in place. An export is what a\n` +
                `world is reconstructed from and nothing keeps a copy of one, so the\n` +
                `terrain would change under everything already derived from it and the\n` +
                `previous bytes would be gone.\n\n` +
                `Export to a NEW directory. If this one holds a failed export, remove it\n` +
                `yourself first -- deliberately, and having looked at what is in it.`);
            process.exit(2);
        }
    }

    const log = args.quiet ? () => {} : (...m) => console.log(...m);

    // The simulation modules log diagnostics to console.log unconditionally.
    // --quiet is meant for scripting, so silence them; warnings and errors
    // still go to stderr.
    const stdout = console.log.bind(console);   // survives the --quiet stub below
    if (args.quiet) console.log = () => {};

    // Imported here rather than at the top so the simulation modules' import-time
    // logging lands after the --quiet stub is installed.
    const { runGeneratePipeline } = await import('../js/pipeline.js');

    const pl = makePlanet(args.planet);
    log(`Generating planet ${seed} — ${args.regions} regions, ${args.plates} plates`);
    if (pl.radiusKm !== 6371 || pl.gravityMS2 !== 9.80665) {
        log(`  ${pl.name}: R=${pl.radiusKm} km (${pl.radiusEarth.toFixed(2)} Earth), `
          + `g=${pl.gravityMS2.toFixed(2)} m/s² (${pl.gravityEarth.toFixed(2)} Earth), `
          + `relief x${pl.reliefScale.toFixed(2)}`);
    }
    const t0 = performance.now();

    let lastPct = -1;
    const ctx = runGeneratePipeline({
        N: args.regions,
        P: args.plates,
        jitter: args.jitter,
        nMag: args.noise,
        numContinents: args.continents,
        continentSizeVariety: args.continentSizeVariety,
        landCoverage: args.landCoverage,
        terrainWarp: args.terrainWarp,
        smoothing: args.smoothing,
        hydraulicErosion: args.hydraulic,
        thermalErosion: args.thermal,
        ridgeSharpening: args.ridge,
        glacialErosion: args.glacial,
        temperatureOffset: args.temperatureOffset,
        precipitationOffset: args.precipitationOffset,
        seed,
        toggledIndices: args.toggledIndices || [],
        skipClimate: !args.climate,
        preserveBasins: args.preserveBasins,
        preserveBasin: args.preserveBasin,
        basinMinDepthKm: args.basinMinDepthKm,
        basinMinAreaKm2: args.basinMinAreaKm2,
        basinMinCells: args.basinMinCells,
        inlandWaterLevels: args.waterLevels,
        preserveBasinList: args.preserveBasinList,
        carveBasinList: args.carveBasinList,
        lithology: args.lithology,
        lithologyStrength: args.lithologyStrength,
        planet: args.planet,
        iceMask: args.iceMask,
    }, (pct, label) => {
        if (pct !== lastPct) { log(`  [${String(pct).padStart(3)}%] ${label}`); lastPct = pct; }
    });

    log(`Generated in ${((performance.now() - t0) / 1000).toFixed(1)}s`);

    if (ctx.glacialSummary && ctx.glacialSummary.glacIdxSource) {
        const g = ctx.glacialSummary;
        log(g.glacIdxSource === 'mask'
            ? `Ice placed from the supplied mask: ${g.glaciatedCells} glaciated regions`
            : `Ice placed by the Earth-calibrated latitude ramp (no --ice-mask): `
              + `${g.glaciatedCells} glaciated regions`);
    }

    if (ctx.hydrology) {
        const big = (ctx.hydrology.dischargeDistribution || [])
            .find(d => d.aboveKm2 === 1e4);
        log(`Hydrology: ${ctx.hydrology.mouthCount} coastal outlets`
          + (big ? ` (${big.outlets} above 10,000 km2)` : '')
          + `, ${ctx.hydrology.oceanBasins.length} marginal seas`);
    }
    if (ctx.lithology) {
        const top = ctx.lithology.composition.slice(0, 3)
            .map(c => `${c.code} ${(100 * c.fraction).toFixed(0)}%`).join(', ');
        log(`Lithology: ${ctx.lithology.exhumedCells} cells exhumed to basement; land surface ${top}`);
    }
    if (ctx.basins) {
        const b = ctx.basins;
        log(`Basins: ${b.catalogue.length} depressions detected, ${b.selected.length} preserved`
          + (b.source && b.source !== 'threshold' ? ` (${b.source})` : ''));
        for (const w of b.warnings) console.warn(`  warning: ${w}`);
    }

    if (args.listBasins) {
        if (!ctx.basins) { console.error('Basin preservation is disabled (--no-basins) — nothing to list.'); return; }
        printBasins(ctx.basins, stdout);
        return;
    }

    // ONE generation, every target. buildExportBundle only READS ctx: it builds
    // a fresh manifest and a fresh file list per call and writes nothing back,
    // so gridding cannot feed into the terrain and every target of one
    // invocation carries the same manifest.hashes.finalElevation by
    // construction rather than by luck. The bundle for one target is released
    // before the next is built, so the peak is the mesh plus one target.
    const multi = args.targets.length > 1;
    for (const target of args.targets) {
        const outDir = target.out ?? defaultOut;
        if (multi) log(`\n── ${gridLabel(target)} → ${outDir}`);

        // buildExportBundle reads the same shape state.curData has, plus the extras
        // the pipeline keeps: plateTable, superPlateData, tectonics.
        const { files, manifest } = buildExportBundle(ctx, {
            raw: target.raw,
            grid: target.grid,
            gridWidth: target.gridWidth,
            gridHeight: target.gridHeight,
            gridMethod: target.gridMethod,
            gridType: target.gridType,
            gridTruncation: target.gridTruncation,
            subgridOrography: target.subgrid,
            groups: target.only,
            onProgress: (frac, label) => log(`  [${String(Math.round(frac * 100)).padStart(3)}%] ${label}`),
        });

        fs.mkdirSync(outDir, { recursive: true });
        let total = 0;

        if (target.zip) {
            const zip = zipStore(files);
            fs.writeFileSync(path.join(outDir, 'planet.zip'), zip);
            total = zip.length;
            log(`\nWrote ${path.join(outDir, 'planet.zip')} — ${humanBytes(total)}, ${files.length} entries`);
        } else {
            for (const f of files) {
                const dest = path.join(outDir, f.name);
                fs.mkdirSync(path.dirname(dest), { recursive: true });
                fs.writeFileSync(dest, f.bytes);
                total += f.bytes.length;
            }
            log(`\nWrote ${files.length} files to ${outDir} — ${humanBytes(total)}`);
        }

        if (target.netcdf) {
            if (!target.grid) {
                console.warn('  --netcdf ignored: NetCDF output is gridded and --no-grid was set');
            } else {
                const n = writeNetCDFFile(outDir, manifest, target, files);
                log(`Wrote ${path.join(outDir, 'planet.nc')} — ${humanBytes(n)}`);
            }
        }

        const nFields = (manifest.grid ?? manifest.raw).fields.length;
        log(`${nFields} fields exported. Field catalogue and units are in manifest.json.`);
    }
}

main().catch(err => {
    console.error(`\nerror: ${err.message}`);
    if (process.env.DEBUG) console.error(err.stack);
    process.exit(1);
});
