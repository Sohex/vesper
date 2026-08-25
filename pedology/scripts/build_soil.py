"""Weather the lithology into soil, under this world's climate.

Between the climate model and the vegetation model there is a step neither of
them does: turning rock into something a plant can root in. Parent material sets
what minerals are available; climate sets how far they have been converted;
relief and erosion set how much regolith survives; and the biosphere sets the
organic fraction, which is what makes this a loop rather than a stage.

Emits the soil map LPJ-GUESS's `SoilInput` reads, one row per land cell:

    Lon Lat sand clay silt orgc ph bulkdensity cn soilc depth awc bedrockfrac andic pfixation

(`SoilInput` skips the last five by name.)

Nothing here is specific to Vesper. Every Earth calibration lives in
`../config/pedogenesis.yaml`; this reads the shared `config/planet.yaml`, the
World Orogen export and an ExoPlaSim climatology, none of which are
world-specific in structure. Porting to another planet is a config change.

    python pedology/scripts/build_soil.py
    python pedology/scripts/build_soil.py --soil-carbon <lpj cpool.out>

The second form closes the loop: on iteration 0 there is no biosphere and soil
carbon is whatever `organic.initial_soil_carbon_kg_m2` declares, and on every
iteration after it is LPJ-GUESS's own answer fed back in.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import netCDF4 as nc
import numpy as np
import yaml

from _paths import (ANALYSIS, CONFIG, COMPONENT_ROOT, DATA, PEDOGENESIS,
                    PROJECT_ROOT, climatology_path)

import climatology as climatology_lib  # noqa: E402  from lib/, via _paths.
# Aliased because `climatology` is a local Path in main(); see build_surface_classes.py.

import builds
from paths import rel  # noqa: E402
import orbit
from gridding import land_fraction_of_class
from orogen import Export, LAND

# Coordinate precision shared with biosphere/scripts/build_lpj_driver.py.
# LPJ-GUESS keys its soil map on an exactly-compared pair of doubles, so both
# files must round identically or every lookup misses.
COORD_DECIMALS = 4

EARTH_YEAR_DAYS = orbit.EARTH_CALENDAR_YEAR_DAYS
KELVIN = 273.15

LAND_COLUMN_CONTRACT = COMPONENT_ROOT / "config" / "land_column_properties.yaml"


def column_base_m() -> float:
    """How deep the simulated land column goes, from the contract that declares it.

    Read rather than written as a literal, because the depth is one fact and
    `pedology/config/land_column_properties.yaml` is where it is declared: the
    physical column is `physical_layer_count` layers of
    `physical_layer_thickness_m`, and the contract's own check refuses a
    declaration whose layers do not reach `column_base_m`.
    """
    decl = yaml.safe_load(LAND_COLUMN_CONTRACT.read_text(encoding="utf-8"))
    return float(decl["geometry"]["column_base_m"])

# Which pH parent group each Orogen rock category maps to. Kept in code rather
# than config because it is a classification of the *export's* vocabulary, not a
# tunable: if Orogen adds a class this must fail loudly, which it does.
PH_GROUP = {
    "morb": "igneous_mafic", "oib": "igneous_mafic",
    "flood_basalt": "igneous_mafic", "arc_basalt": "igneous_mafic",
    "arc_andesite": "igneous_felsic", "rift_bimodal": "igneous_felsic",
    "granite": "igneous_felsic", "granodiorite": "igneous_felsic",
    "gneiss": "metamorphic", "schist": "metamorphic",
    "quartzite": "igneous_felsic", "melange": "metamorphic",
    "shelf_clastic": "sedimentary_clastic",
    "foreland_clastic": "sedimentary_clastic",
    "continental_clastic": "sedimentary_clastic",
    "pelagic": "sedimentary_clastic",
    "carbonate": "carbonate",
    "evaporite": "evaporite", "playa_clastic": "evaporite",
    "water": "sedimentary_clastic",
}


def lithology_fractions(config: dict) -> tuple[dict[str, np.ndarray], Export, Path]:
    """Area share of each rock class within every cell's land, from the mesh.

    Fractions rather than the dominant class, because texture is continuous and
    a cell that is half granite and half basalt really does carry a mixture.
    That is the opposite of the soil-code case, where averaging categories would
    have been meaningless.
    """
    mesh = Export(builds.mesh_export(config))
    resolution = str(config["model"]["resolution"]).upper()
    grid_dir = builds.grid_export(config, resolution)
    rock = mesh.substrate_class
    classes = {c["id"]: c["code"] for c in mesh.manifest["lithology"]["rockClasses"]}

    fractions: dict[str, np.ndarray] = {}
    for rock_id, code in classes.items():
        selected = rock == rock_id
        if selected.any():
            fractions[code] = land_fraction_of_class(mesh, grid_dir, selected)
    return fractions, mesh, grid_dir


def subgrid_slope(mesh: Export, grid_dir: Path, baseline_km: float) -> np.ndarray:
    """Within-cell elevation spread, expressed over a DECLARED run.

    Area-weighted standard deviation of mesh-region elevation inside each grid
    cell, over `baseline_km` of ground. Land regions only, so a coastal cell is
    not handed the continental shelf as a hillslope.

    THE RUN IS DECLARED AND IS NOT THE MESH SPACING. Dividing by the export's
    own `avgEdgeKm` put the region count into the answer: the divisor halves
    when the region count quadruples, so the quotient doubles on terrain that
    has converged. Measured between this project's two builds, WORLD-UYFB and
    `analysis/subgrid_slope_support.py`: the quotient moves by 2.01 to 2.14x
    across the land quantiles while the spread it is built from moves by 1.005
    to 1.069x, and the residual once the spacing ratio is divided out is 1.000
    at every quantile. The length in the definition was the whole of the shift.
    That is the compound topographic index's mechanism, not the self-affine one
    the scarp gradient has and not a population bias; over a declared run the
    same quantiles agree to 1.069x, inside the 1.15x bar
    `notes/audits/orogen-resolution.md` fixed for this class. The baseline is a
    unit for the spread and not a length this statistic samples at, so it does
    not have to be the 90 km the scarp relief term declares.

    WHAT THE MAGNITUDE IS. A spread over a declared length, and not the
    terrain's hillslope gradient: a T42 cell is hundreds of kilometres across
    and real catenas run at 100 m. The field carries the PATTERN -- mountains
    thin, basins deep -- and `catena.slope_transport` is declared against real
    hillslope gradients rather than fitted to this distribution. So no absolute
    threshold keyed to a measured gradient belongs on it, and what survives a
    change of baseline is the ordering of the cells rather than the value.
    """
    from gridding import cell_moments, region_cells

    cell, nlat, nlon = region_cells(mesh, grid_dir)
    elevation_m = mesh.elevation_km.astype(np.float64) * 1000.0
    is_land = mesh.surface_class == LAND
    _, variance, _, _ = cell_moments(cell, nlat * nlon,
                                     mesh.cell_area.astype(np.float64),
                                     elevation_m, population=is_land)
    return (np.sqrt(variance) / (baseline_km * 1000.0)).reshape(nlat, nlon)


def weathering_intensity(runoff_mm_yr: np.ndarray, temperature_c: np.ndarray,
                         params: dict, reference: float | None = None) -> np.ndarray:
    """Walker-Hays-Kasting weathering intensity, normalised to Earth land means.

    Dimensionless and clipped. It is an intensity rather than a rate: the time
    integral is folded into the reference, because soil age is not known on
    either world well enough to carry explicitly.
    """
    if reference is None:
        reference = params["reference_runoff_mm_per_earth_year"]
    q_ratio = np.maximum(runoff_mm_yr, 0.0) / reference
    thermal = np.exp((temperature_c - params["reference_temperature_c"])
                     / params["temperature_e_folding_k"])
    intensity = np.power(np.maximum(q_ratio, 1e-6), params["runoff_exponent"]) * thermal
    return np.clip(intensity, params["minimum"], params["maximum"])


def weather_texture(fractions: dict[str, np.ndarray], intensity: np.ndarray,
                    texture_cfg: dict) -> dict[str, np.ndarray]:
    """Mix parent textures, then convert weatherable minerals to clay.

    Quartz is inert and never becomes clay, which is the whole reason granite and
    basalt diverge under the same climate: both start sandy, but only one of them
    has sand that can weather.
    """
    shape = intensity.shape
    sand = np.zeros(shape)
    silt = np.zeros(shape)
    clay = np.zeros(shape)
    quartz = np.zeros(shape)
    total = np.zeros(shape)

    for code, share in fractions.items():
        entry = texture_cfg.get(code)
        if entry is None:
            raise SystemExit(
                f"rock class {code!r} has no texture in pedogenesis.yaml. Add it "
                f"and record the reasoning.")
        sand += share * entry["sand"]
        silt += share * entry["silt"]
        clay += share * entry["clay"]
        quartz += share * entry["quartz"]
        total += share

    # Cells with no mesh land carry no lithology; leave them at the land mean
    # later rather than dividing by zero here.
    covered = total > 1e-9
    for array in (sand, silt, clay, quartz):
        array[covered] /= total[covered]

    # Only part of the weatherable fraction becomes clay-sized silicate; the rest
    # leaves as solute or becomes sesquioxide, neither of which a texture
    # analysis counts. Calibrated against SoilGrids over sixteen type
    # localities; see pedogenesis.yaml.
    weatherable = np.clip(1.0 - quartz - clay, 0.0, 1.0) * texture_cfg["clay_yield"]
    converted = weatherable * (1.0 - np.exp(-texture_cfg["clay_conversion"] * intensity))

    ratio = texture_cfg["sand_to_silt_loss_ratio"]
    denominator = np.maximum(sand * ratio + silt, 1e-9)
    from_sand = np.minimum(converted * sand * ratio / denominator, sand)
    from_silt = np.minimum(converted * silt / denominator, silt)

    return {
        "sand": np.clip(sand - from_sand, 0.0, 1.0),
        "silt": np.clip(silt - from_silt, 0.0, 1.0),
        "clay": np.clip(clay + from_sand + from_silt, 0.0, 1.0),
        "quartz": quartz,
        "covered": covered,
    }


def andisol_properties(fractions: dict[str, np.ndarray], intensity: np.ndarray,
                       precipitation_mm_yr: np.ndarray, cfg: dict
                       ) -> dict[str, np.ndarray]:
    """Andic and vitric fractions, and what they do to the soil.

    Two gates multiply. `resupply` is the area of the cell receiving ongoing
    volcanic ejecta -- the arc, and only the arc, because it is placed as a band
    about the volcanic front and is the one part of this pipeline that knows
    volcanism is still happening. `allophanic` is whether leaching is strong
    enough that the glass weathers to allophane rather than to halloysite, at
    Dahlgren's ~1500 mm precipitation boundary.

    Development then splits the resupplied area between vitric (glass still
    dominant, weak andic expression) and andic (allophane formed). Halloysitic
    volcanic terrain is neither: it is an ordinary clay soil, and on this world
    it is most of the volcanic terrain.
    """
    zero = np.zeros_like(intensity)
    resupply = sum((fractions.get(c, zero) for c in cfg["resupplied_classes"]),
                   start=np.zeros_like(intensity))
    # Soft boundary, because the source says "generally less than about".
    allophanic = 1.0 / (1.0 + np.exp(
        -(precipitation_mm_yr - cfg["allophane_precipitation_mm"])
        / cfg["allophane_transition_mm"]))
    # Vitric -> andic with weathering progress. Saturating, so it approaches
    # full andic expression rather than crossing a step.
    development = intensity / (intensity + cfg["development_half_saturation_w"])
    andic = resupply * allophanic * development
    vitric = resupply * allophanic * (1.0 - development)
    # Volcanic terrain that weathered to halloysite instead. Tracked because it
    # is the answer to "why are there so few andisols here", and because it is
    # where silica stays in the profile rather than leaving in solution.
    halloysitic = resupply * (1.0 - allophanic)
    # The ANDIC contribution to phosphate fixation, not a soil's absolute
    # retention. Halloysitic and non-volcanic material contributes zero here
    # because it has no andic fixation, which is not the same as saying it
    # retains no phosphorus -- it retains it the ordinary way, which every other
    # soil in this map does too and which this term is not about. Read it as:
    # the fraction of released P in this cell that andic material takes out of
    # circulation, over and above whatever a normal soil would do.
    fixation = (cfg["phosphate_retention_andic"] * andic
                + cfg["phosphate_retention_vitric"] * vitric)
    return {"andic": andic, "vitric": vitric, "halloysitic": halloysitic,
            "resupplied": resupply, "andic_p_fixation": fixation}


def regolith_depth(intensity: np.ndarray, relief_m: np.ndarray,
                   runoff_mm_yr: np.ndarray, erodibility: np.ndarray,
                   params: dict, weathering_ref: float) -> np.ndarray:
    """Steady-state regolith thickness, production against erosion, bounded.

    Heimsath's exponential soil production function gives a steady state of
    `h = h_star * ln(P0 / E)`, and that is what this used to compute. It has a
    disqualifying property for a whole planet: it diverges as erosion approaches zero and
    goes to minus infinity as erosion grows, so both ends had to be clipped, and
    a quarter of this world's land sat on each clip. That is not a bimodal
    planet, it is a model railing.

    The reason is dynamic range. Filling 0 to 5 m through a logarithm with
    `h_star` at Heimsath's 0.5 m needs production over erosion to span e^10,
    about 22,000. Nothing in these inputs spans that, so almost every cell lands
    outside and gets clipped.

    Replaced by a saturating form:

        depth = maximum_depth * P / (P + erosion_weight * E)

    Bounded at both ends by construction and monotone in the right directions.
    As erosion vanishes the profile approaches `maximum_depth`, which is the
    physical statement that a weathering front cannot advance forever because
    water and oxygen have to reach it through what has already accumulated. As
    erosion grows the profile thins smoothly to bare rock.

    This is a parameterisation and not a derivation, unlike the Heimsath form it
    replaces, and that is the trade: an honest curve that spans the range against
    a principled one that cannot be evaluated over it.
    """
    # Moisture-driven denudation plus a floor that does not need runoff. A dry
    # slope still loses material to wind, dry ravel and creep; without the floor
    # erosion is exactly zero wherever P - E is, and every arid cell pins to the
    # ceiling.
    moisture = (np.maximum(runoff_mm_yr, 0.0) / weathering_ref
                + params["dry_erosion_baseline"])
    erosion = (erodibility
               * np.maximum(relief_m, 0.0) / params["erosion_reference_relief_m"]
               * moisture)
    production = np.maximum(intensity, 1e-9)
    depth = (params["maximum_depth_m"] * production
             / (production + params["erosion_weight"] * erosion))
    return np.clip(depth, params["minimum_depth_m"], params["maximum_depth_m"])


def regolith_depth_expectation(mesh: Export, grid_dir: Path,
                               intensity: np.ndarray, relief_m: np.ndarray,
                               runoff_mm_yr: np.ndarray, params: dict,
                               weathering_ref: float):
    """The depth law evaluated on EACH ROCK, then area-averaged over the cell.

    The law is convex in the erodibility, which spans a factor of fourteen
    across the export's lithology table, so mixing a cell's rocks into one
    erodibility and then applying the law is not the same operation as applying
    the law to each rock and averaging the depths. It is smaller, always:
    `notes/audits/nonlinear-spatial-reductions.md` section 1 sizes the
    difference at 0.19 to 0.28 m of regolith at the land mean against the 0.02 m
    `regolith.minimum_depth_m` bar, one-signed over the whole sweep of the law's
    one free ratio and at every rung of the ladder, and refinement shrinks it
    without removing it.

    The reduction is `gridding.cell_expectation` and is CALLED rather than
    written again here. That operator takes the LAW as a callable, so a caller
    who has already reduced a field cannot hand the reduced field in by
    accident; the order is a property of the call site and not of whoever last
    edited it.

    Erodibility is the only one of the law's inputs the mesh resolves. The
    others are cell fields -- weathering intensity and runoff come from the
    climatology's own grid and relief is a difference between neighbouring cell
    centres -- so each region sees its own rock under its cell's climate, which
    is `hydrography/notes/subgrid-water-table.md` section 2: sub-grid
    information reaches a cell-scale parameter as a statistic of the cell's own
    distribution and never as a resolved position.

    Returns `(depth, depth_at_the_mixed_erodibility, erodibility, ledger)`,
    all four on the grid except the ledger. The second is what the mixed order
    would have given and is kept so the report can carry the gap rather than
    leave a later reader to re-derive it.
    """
    from gridding import cell_expectation, region_cells, transfer_ledger

    cell, nlat, nlon = region_cells(mesh, grid_dir)
    ncell = nlat * nlon
    if intensity.shape != (nlat, nlon):
        raise SystemExit(
            f"the climatology is {intensity.shape} and the grid export is "
            f"{(nlat, nlon)}. Both must be the configured rung; see lib/rungs.py.")
    area = mesh.cell_area.astype(np.float64)
    is_land = mesh.surface_class == LAND

    entries = mesh.manifest["lithology"]["rockClasses"]
    erod_of_class = np.ones(max(int(e["id"]) for e in entries) + 1)
    for entry in entries:
        erod_of_class[int(entry["id"])] = float(entry["erodibility"])
    erod_region = erod_of_class[mesh.substrate_class.astype(int)]

    def at_region(field: np.ndarray) -> np.ndarray:
        """A cell field, as every region inside that cell sees it."""
        return np.asarray(field, dtype=np.float64).reshape(-1)[cell]

    region_intensity = at_region(intensity)
    region_relief = at_region(relief_m)
    region_runoff = at_region(runoff_mm_yr)

    def law(erodibility: np.ndarray) -> np.ndarray:
        return regolith_depth(region_intensity, region_relief, region_runoff,
                              erodibility, params, weathering_ref)

    depth, mixed, covered = cell_expectation(cell, ncell, area, law,
                                             erod_region, is_land)
    # A cell holding no land region of the mesh has no rocks to average over.
    # It keeps the erodibility of 1.0 this file has always given it, and the
    # law at that erodibility is then its depth, so the substitution is visible
    # as a value rather than as a zero from an empty bin.
    mixed = np.where(covered, mixed, 1.0).reshape(nlat, nlon)
    at_mixed = regolith_depth(intensity, relief_m, runoff_mm_yr, mixed,
                              params, weathering_ref)
    depth = np.where(covered.reshape(nlat, nlon), depth.reshape(nlat, nlon),
                     at_mixed)
    return depth, at_mixed, mixed, transfer_ledger(cell, ncell, area, is_land)


def soil_ph(fractions: dict[str, np.ndarray], runoff_mm_yr: np.ndarray,
            endorheic: np.ndarray, params: dict, reference_runoff: float
            ) -> np.ndarray:
    """Parent pH, leached down by drainage, pushed up where drainage is closed.

    Deliberately takes runoff, not whatever `weathering.moisture_variable`
    selects. Leaching is base cations physically leaving the profile, which
    requires water to drain through it; rain that falls and evaporates carries
    nothing away. So this stays on runoff even when weathering is driven by
    precipitation, and pH is consequently the one soil property that does not
    move when that switch is flipped. Intended, not an oversight.
    """
    shape = runoff_mm_yr.shape
    parent = np.zeros(shape)
    total = np.zeros(shape)
    for code, share in fractions.items():
        group = PH_GROUP.get(code)
        if group is None:
            raise SystemExit(f"rock class {code!r} has no pH group; add it to PH_GROUP")
        parent += share * params["parent_by_category"][group]
        total += share
    covered = total > 1e-9
    parent[covered] /= total[covered]
    parent[~covered] = params["parent_by_category"]["sedimentary_clastic"]

    leaching = np.log1p(np.maximum(runoff_mm_yr, 0.0) / reference_runoff)
    ph = parent - params["leaching_slope"] * leaching
    ph += endorheic * params["endorheic_alkalinity_bonus"]
    return np.clip(ph, params["minimum"], params["maximum"])


def organic_properties(soil_carbon_kg_m2: np.ndarray, params: dict
                       ) -> tuple[np.ndarray, np.ndarray]:
    """Organic mass fraction and bulk density from soil carbon.

    Bulk density uses the standard reciprocal mixing rule (Adams 1973), which is
    why organic soils are light: the two components occupy volume additively, not
    mass additively.
    """
    organic_matter = (np.maximum(soil_carbon_kg_m2, 0.0)
                      / params["carbon_fraction_of_organic_matter"])
    # Solve the mixing rule for the mass fraction that this areal density implies
    # at the reference depth, iterating once from a mineral-only first guess.
    depth = params["reference_depth_m"]
    bulk = np.full(soil_carbon_kg_m2.shape, float(params["mineral_bulk_density_kg_m3"]))
    fraction = np.zeros_like(bulk)
    for _ in range(8):
        fraction = np.clip(organic_matter / np.maximum(bulk * depth, 1e-9),
                           0.0, params["maximum_organic_fraction"])
        inverse = (fraction / params["organic_bulk_density_kg_m3"]
                   + (1.0 - fraction) / params["mineral_bulk_density_kg_m3"])
        bulk = 1.0 / np.maximum(inverse, 1e-12)
    return fraction, bulk


def read_soil_carbon(path: Path, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    """Read LPJ-GUESS soil carbon (cpool.out) onto the grid, last year per cell."""
    lines = path.read_text().splitlines()
    header = lines[0].split()
    try:
        column = header.index("SoilC")
    except ValueError:
        raise SystemExit(f"{path} has no SoilC column; is it a cpool.out?")
    # LPJ-GUESS prints coordinates at two decimals in its output files, which is
    # coarser than the driver's four. At any sane model resolution two decimals
    # is still unique per cell, so match on that rather than on the grid's own
    # precision, which would never agree.
    output_decimals = 2

    latest: dict[tuple[float, float], tuple[int, float]] = {}
    for line in lines[1:]:
        parts = line.split()
        if len(parts) <= column:
            continue
        key = (round(float(parts[0]), output_decimals),
               round(float(parts[1]), output_decimals))
        year = int(float(parts[2]))
        if key not in latest or year > latest[key][0]:
            latest[key] = (year, float(parts[column]))
    grid = np.zeros((len(lat), len(lon)))
    lon_signed = np.where(lon > 180.0, lon - 360.0, lon)
    hits = 0
    for j, cell_lat in enumerate(lat):
        for i, cell_lon in enumerate(lon_signed):
            key = (round(float(cell_lon), output_decimals),
                   round(float(cell_lat), output_decimals))
            if key in latest:
                grid[j, i] = latest[key][1]
                hits += 1
    if hits == 0:
        raise SystemExit(
            f"{path} shares no coordinates with the grid; are they the same run?")
    print(f"soil carbon: matched {hits} cells from {path.name}")
    return grid


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--climatology", type=Path, default=None)
    parser.add_argument("--soil-carbon", type=Path, default=None,
                        help="LPJ-GUESS cpool.out from the previous iteration. "
                             "Omit for iteration 0, which has no biosphere.")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--iteration", type=int, default=None,
                        help="loop iteration number, for the provenance record")
    args = parser.parse_args()

    config = yaml.safe_load(CONFIG.read_text())
    pedo = yaml.safe_load(PEDOGENESIS.read_text())
    # Resolve before use, not just before storing. Scripts here anchor their
    # paths from the file location rather than the cwd, so a `--climatology`
    # given relative to wherever the caller stood has to be made absolute before
    # anything opens it.
    climatology = (args.climatology or climatology_path()).resolve()
    if not climatology.is_file():
        raise SystemExit(f"{climatology} does not exist")
    # Deliberate, not assumed: this soil map pairs a climatology with a
    # lithology, and they have to be the same world.
    from provenance import require_build
    require_build(climatology, "climatology", config)

    with nc.Dataset(climatology) as data:
        lat = np.asarray(data["lat"][:], dtype=float)
        lon = np.asarray(data["lon"][:], dtype=float)
        land = np.asarray(data["lsm"][0], dtype=float) > 0.5
        temperature = climatology_lib.annual_mean_of(data, "tas") - KELVIN
        # Runoff and precipitation are rates; annualise on Earth years so the
        # Earth-calibrated weathering law is fed the units it was fitted in.
        # Runoff, from whichever source the config trusts. See pedogenesis.yaml:
        # the reported mrro field and the land water budget disagree by 6.6x and
        # the budget wins.
        scale = 1000.0 * 86400.0 * EARTH_YEAR_DAYS
        source = pedo["weathering"].get("runoff_source", "p_minus_e")
        if source == "mrro":
            runoff = climatology_lib.annual_mean_of(data, "mrro") * scale
        elif source == "p_minus_e":
            evaporation = -climatology_lib.annual_mean_of(data, "evap")
            runoff = (climatology_lib.annual_mean_of(data, "pr")
                      - evaporation) * scale
            # Clipped at zero. Per cell the residual can go negative where the
            # postprocessing is noisy or a river routes water through; at the
            # land mean it is the conserved quantity and no clipping applies.
            runoff = np.maximum(runoff, 0.0)
        else:
            raise SystemExit(
                f"weathering.runoff_source is {source!r}; expected 'p_minus_e' "
                f"or 'mrro'")
        precip = climatology_lib.annual_mean_of(data, "pr") \
            * 1000.0 * 86400.0 * EARTH_YEAR_DAYS
        elevation = climatology_lib.annual_mean_of(data, "sg") \
            / float(config["planet"]["gravity_m_s2"])
        # Timestep extrema, so these carry the full 30-hour diurnal swing and a
        # bin can straddle freezing even where the bin mean never does.
        bin_min = np.asarray(data["mint"][:], dtype=float) - KELVIN
        bin_max = np.asarray(data["maxt"][:], dtype=float) - KELVIN

    fractions, mesh, grid_dir = lithology_fractions(config)

    # Local relief, as the spread of surface height across each cell's neighbours.
    # A stand-in for slope that needs no extra field and no mesh gradient.
    relief = np.zeros_like(elevation)
    for axis in (0, 1):
        forward = np.roll(elevation, 1, axis=axis)
        backward = np.roll(elevation, -1, axis=axis)
        relief = np.maximum(relief, np.maximum(np.abs(elevation - forward),
                                               np.abs(elevation - backward)))

    # Slope for the transport term, from the native mesh rather than from this
    # grid, and it has to be.
    #
    # Catena is a hillslope process and a T42 cell is hundreds of kilometres
    # across. Differencing neighbouring cell centres gives a land-mean gradient
    # three orders of magnitude below a real hillslope, and the term does
    # nothing. The mesh resolves the within-cell elevation spread, which is a
    # far better representative statistic; the run it is quoted over comes from
    # `catena.gradient_baseline_km` and NOT from the mesh, for the reason
    # `subgrid_slope` gives.
    #
    # It is still an underestimate: real catenas run at 100 m scale and gradients
    # of 0.1 to 0.5. So this term reproduces the *pattern*, mountains thin and
    # basins deep, rather than the absolute magnitude, and `slope_transport` is
    # declared against real hillslope gradients. Recorded rather than tuned away.
    tan_beta = subgrid_slope(mesh, grid_dir,
                             float(pedo["catena"]["gradient_baseline_km"]))

    # Fraction of the orbit whose diurnal range straddles freezing. Freeze-thaw
    # cycling shatters rock; ground frozen solid all bin does not.
    frost_fraction = np.mean((bin_min < 0.0) & (bin_max > 0.0), axis=0)

    endorheic = land_fraction_of_class(mesh, grid_dir, mesh.is_endorheic.astype(bool))

    # Both moisture drivers, so the bracket is visible whichever one is selected.
    intensity_by_runoff = weathering_intensity(runoff, temperature, pedo["weathering"])
    # Against Earth's precipitation, not Earth's runoff. The two branches have to
    # be normalised against their own reference or the comparison is a unit error
    # rather than a bracket.
    intensity_by_precip = weathering_intensity(
        precip, temperature, pedo["weathering"],
        pedo["weathering"]["reference_precipitation_mm_per_earth_year"])
    selector = pedo["weathering"].get("moisture_variable", "runoff")
    if selector == "runoff":
        intensity = intensity_by_runoff
    elif selector == "precipitation":
        intensity = intensity_by_precip
    else:
        raise SystemExit(
            f"weathering.moisture_variable is {selector!r}; expected "
            f"'runoff' or 'precipitation'")
    texture = weather_texture(fractions, intensity, pedo["texture"])
    # Erodibility is NOT mixed before the law is applied; see
    # `regolith_depth_expectation` for the operator and for what the two orders
    # are worth on this map.
    depth, depth_at_mixed, erodibility, regolith_ledger = \
        regolith_depth_expectation(
            mesh, grid_dir, intensity, relief, runoff, pedo["regolith"],
            pedo["weathering"]["reference_runoff_mm_per_earth_year"])

    # Catena. Frost shattering adds production; slope takes it away again. The
    # chemical side is untouched: clay fraction still comes from the weathering
    # intensity, and this only decides how much material there is and how much of
    # the fine fraction stayed put.
    catena = pedo["catena"]
    # What the mixed order would have cost, in the law's own units, before the
    # catena factors: they are cell-scale statistics and multiply both orders by
    # the same number, so the whole of the difference is the operator's.
    regolith_gap = depth - depth_at_mixed
    # A PRE-REGISTERED INVARIANT, not a comparison. The depth law is a single
    # saturating function of the erodibility and its curvature does not change
    # sign, so the expectation cannot come out below the law at the mean: the
    # mixed order can only ever make the profile thinner. A negative gap past
    # round-off therefore says the law, the population or the binning moved,
    # and it is a failure with a right answer rather than a difference to note.
    if float(regolith_gap.min()) < -1e-9:
        raise SystemExit(
            f"the regolith expectation came out {float(regolith_gap.min()):.6g} m "
            "BELOW the law at the mixed erodibility. The depth law is convex in "
            "the erodibility, so that cannot happen while the two arms share a "
            "law, a population and a binning; one of the three has moved.")
    catena_factor = ((1.0 + catena["frost_production_bonus"] * frost_fraction)
                     / (1.0 + catena["slope_transport"] * tan_beta))
    depth = depth * catena_factor
    depth_at_mixed = depth_at_mixed * catena_factor
    # The one place after the law where the two orders could part again, since
    # a clip is not linear. Counted rather than argued about.
    catena_clipped = ((depth < pedo["regolith"]["minimum_depth_m"])
                      | (depth > pedo["regolith"]["maximum_depth_m"]))
    depth = np.clip(depth, pedo["regolith"]["minimum_depth_m"],
                    pedo["regolith"]["maximum_depth_m"])
    depth_at_mixed = np.clip(depth_at_mixed, pedo["regolith"]["minimum_depth_m"],
                             pedo["regolith"]["maximum_depth_m"])

    # Fines are shed preferentially downhill, so a slope keeps the coarse
    # fraction. Moves clay to sand, conserving the total.
    fines_loss = np.clip(catena["slope_fines_loss"] * tan_beta,
                         0.0, catena["maximum_fines_loss"])
    moved = texture["clay"] * fines_loss
    texture["clay"] = texture["clay"] - moved
    texture["sand"] = texture["sand"] + moved
    ph = soil_ph(fractions, runoff, endorheic, pedo["ph"],
                 pedo["weathering"]["reference_runoff_mm_per_earth_year"])

    # Andisols. Volcanism as a process rather than a composition: see the
    # `andisol` block in pedogenesis.yaml for why this is gated on the arc alone
    # and on a precipitation threshold rather than on rock type.
    andisol = andisol_properties(fractions, intensity, precip, pedo["andisol"])

    if args.soil_carbon:
        carbon = read_soil_carbon(args.soil_carbon, lon, lat)
    else:
        carbon = np.full(temperature.shape,
                         float(pedo["organic"]["initial_soil_carbon_kg_m2"]))
    organic_fraction, bulk_density = organic_properties(carbon, pedo["organic"])

    # Andic material is far less dense than its texture and organic content
    # imply, because the porosity is inside allophane and humus rather than
    # between grains. Blended by andic fraction, not substituted, since a cell
    # is only partly arc terrain.
    bulk_density = (bulk_density * (1.0 - andisol["andic"])
                    + pedo["andisol"]["bulk_density_andic"] * andisol["andic"])

    # Plant-available water capacity, mm: volumetric capacity from texture times
    # the depth of regolith that actually exists, cut at the declared base of
    # the simulated land column.
    #
    # THE CUT IS THE DEFINITION AND NOT A TRUNCATION. This is a PLANT-AVAILABLE
    # capacity. The land column property contract declares fifteen 100 mm
    # layers to `column_base_m`, LPJ-GUESS's root distribution spans exactly
    # those fifteen layers, and the contract's rootable base is shallower still
    # at `rootable_base_m`, so no root and no evaporating surface in this
    # pipeline can reach water below the column base. Regolith runs deeper than
    # that on part of this map, and multiplying a plant-available capacity by
    # the whole of it puts water that nothing can draw on into a number that
    # says how much can be drawn on.
    #
    # The water below the base is not lost physics. It belongs to the aquifer
    # term, which the contract already declares this column stops at, and
    # `hydrography/config/land_water_ledger.yaml` registers
    # `transient_saturated_storage` as an UNREPRESENTABLE owned by PLHY-4
    # because the groundwater solver is steady-state and carries no storage
    # change. Deepening the column instead -- 26 layers to cover the regolith
    # p90 and 43 for the p99, against NSOILLAYER appearing across six
    # LPJ-GUESS modules -- becomes worth doing when that solver gains storage,
    # and not before. WORLD-7702.
    water = pedo["water"]
    column_base = column_base_m()
    column_depth = np.minimum(depth, column_base)
    volumetric = (water["volumetric_capacity_by_texture"]["sand"] * texture["sand"]
                  + water["volumetric_capacity_by_texture"]["silt"] * texture["silt"]
                  + water["volumetric_capacity_by_texture"]["clay"] * texture["clay"]
                  + water["volumetric_capacity_organic"] * organic_fraction
                  # Noncrystalline material holds water the texture terms cannot
                  # see. Additive for the same reason the organic term is.
                  + pedo["andisol"]["volumetric_capacity_allophane"]
                  * andisol["andic"])
    water_capacity = np.clip(volumetric * column_depth * 1000.0,
                             water["minimum_mm"], water["maximum_mm"])

    # The same capacity at both ends of the declared brackets on the volumetric
    # constants. These are DECLARED values rather than measured ones -- Saxton
    # and Rawls (2006) fixes the shape but excludes the pure endmembers this
    # mixes between -- and this field sets ExoPlaSim's dwmax, whose overflow IS
    # its runoff, so reporting one number for it would hide what the declaration
    # is worth. LITH-24, and the same treatment aeolian/config/dust.yaml gets.
    def _capacity_at(end: int) -> np.ndarray:
        bracket = water["volumetric_capacity_by_texture_bracket"]
        vol = (bracket["sand"][end] * texture["sand"]
               + bracket["silt"][end] * texture["silt"]
               + bracket["clay"][end] * texture["clay"]
               + water["volumetric_capacity_organic_bracket"][end] * organic_fraction
               + pedo["andisol"]["volumetric_capacity_allophane"] * andisol["andic"])
        return np.clip(vol * column_depth * 1000.0,
                       water["minimum_mm"], water["maximum_mm"])

    water_capacity_low = _capacity_at(0)
    water_capacity_high = _capacity_at(1)

    # Plant-available water below the bedrock contact, as a fraction of what the
    # soil above holds per unit volume. A function of weathering, because that
    # is what turns impermeable rock into saprock and then saprolite. See the
    # anchors in pedogenesis.yaml; this is not a small correction.
    bedrock = pedo["regolith"]["bedrock_water"]
    bedrock_fraction = (bedrock["minimum"]
                        + (bedrock["maximum"] - bedrock["minimum"])
                        * (1.0 - np.exp(-bedrock["shape"] * intensity)))

    DATA.mkdir(parents=True, exist_ok=True)
    # NOTE: a relative --output resolves against the caller's cwd, not this
    # file's anchor; `rel` records an out-of-tree path as absolute rather than
    # raising.
    # Per build. Soil texture derives from lithology, so a soil map belongs to
    # the terrain it was computed from, and LPJ-GUESS eats this file directly.
    if args.output is None:
        import sys as _sys
        _sys.path.insert(0, str(PROJECT_ROOT / "lib"))
        import builds as _b
        default_out = _b.soilmap(config)
        default_out.parent.mkdir(parents=True, exist_ok=True)
    else:
        default_out = args.output
    output = Path(default_out).resolve()
    lon_signed = np.where(lon > 180.0, lon - 360.0, lon)
    rows = np.argwhere(land)
    with output.open("w") as handle:
        # `depth` and `awc` are extra columns LPJ-GUESS's own SoilInput ignores
        # (it matches columns by name and skips what it does not know). They are
        # here so one artifact carries the whole soil: biosphere reads depth to
        # scale water capacity, and exoplasim reads awc to set its bucket.
        # `andic` and `pfixation` are two more columns LPJ-GUESS's SoilInput
        # skips by name. They are here because the C-N-P fork needs them: andic
        # material fixes phosphorus rather than supplying it, and no other
        # column in this file carries that.
        handle.write("Lon Lat sand clay silt orgc ph bulkdensity cn soilc "
                     "depth awc bedrockfrac andic pfixation\n")
        for j, i in rows:
            handle.write(
                f"{lon_signed[i]:.{COORD_DECIMALS}f} {lat[j]:.{COORD_DECIMALS}f} "
                f"{texture['sand'][j, i]:.4f} {texture['clay'][j, i]:.4f} "
                f"{texture['silt'][j, i]:.4f} {organic_fraction[j, i]:.5f} "
                f"{ph[j, i]:.3f} {bulk_density[j, i]:.1f} "
                f"{pedo['organic']['carbon_nitrogen_ratio']:.1f} "
                f"{carbon[j, i]:.4f} "
                f"{depth[j, i]:.4f} {water_capacity[j, i]:.2f} "
                f"{bedrock_fraction[j, i]:.4f} "
                f"{andisol['andic'][j, i]:.4f} "
                f"{andisol['andic_p_fixation'][j, i]:.4f}\n")

    weights = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, len(lon)))
    lw = weights[land]

    def mean(field: np.ndarray) -> float:
        return float(np.average(field[land], weights=lw))

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "iteration": args.iteration,
        "closes_loop_with": (str(args.soil_carbon) if args.soil_carbon
                             else "nothing; iteration 0 has no biosphere"),
        "climatology": rel(climatology),
        "climatology_sha256": hashlib.sha256(climatology.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
        "pedogenesis_sha256": hashlib.sha256(PEDOGENESIS.read_bytes()).hexdigest(),
        "source_build": config.get("source_build"),
        "terrain_hash": mesh.terrain_hash,
        "orbital_year_earth_days": orbit.orbital_year_days(config),
        "land_cells": int(len(rows)),
        "moisture_variable": selector,
        "runoff_source": pedo["weathering"].get("runoff_source", "p_minus_e"),
        "andisols": {
            "note": ("Andic properties need ONGOING ejecta supply and enough "
                     "leaching to weather glass to allophane rather than "
                     "halloysite. The arc is resupplied because subduction "
                     "is ongoing; flood_basalt and oib are excluded on a "
                     "duration argument -- a province erupts for order 1% "
                     "of its life -- not for want of an age field. "
                     "Fractions are of total land area."),
            "resupplied_land_fraction": mean(andisol["resupplied"]),
            "andic_land_fraction": mean(andisol["andic"]),
            "vitric_land_fraction": mean(andisol["vitric"]),
            "halloysitic_land_fraction": mean(andisol["halloysitic"]),
            "unresupplied_land_fraction": mean(sum(
                (fractions.get(c, np.zeros_like(intensity))
                 for c in pedo["andisol"]["unresupplied_classes"]),
                start=np.zeros_like(intensity))),
            "land_mean_andic_p_fixation": mean(andisol["andic_p_fixation"]),
            "andic_p_fixation_note": (
                "Areal fraction of released P that andic material fixes, over "
                "and above ordinary soil retention. Zero off andic ground by "
                "construction; that is not a claim those soils retain no P."),
            "allophane_precipitation_mm": pedo["andisol"][
                "allophane_precipitation_mm"],
        },
        "weathering_bracket": {
            "note": ("Weathering intensity under each moisture driver, land "
                     "means. The spread is the honest uncertainty in how "
                     "weathered this world's soils are, and it is wide because "
                     "the climate model's land hydrology is a uniform bucket."),
            "by_runoff": mean(intensity_by_runoff),
            "by_precipitation": mean(intensity_by_precip),
            "runoff_ratio": mean(runoff) / max(mean(precip), 1e-9),
            "earth_land_runoff_ratio_for_reference": 0.35,
        },
        "regolith_aggregation": {
            "operator": ("gridding.cell_expectation, NONLINEAR, over the land "
                         "population of the mesh: the depth law is evaluated on "
                         "each region's own rock and the depths are then "
                         "area-averaged. lib/gridding.py owns the reduction "
                         "operators and this file calls one rather than "
                         "reimplementing it."),
            "law": ("depth = maximum_depth * P / (P + erosion_weight * E), with "
                    "E proportional to the region's erodibility"),
            "population": "surface_class == LAND",
            "subgrid_input": "erodibility, from the export's substrate class",
            "bar_m": pedo["regolith"]["minimum_depth_m"],
            "land_mean_erodibility": mean(erodibility),
            "land_mean_depth_at_mixed_erodibility_m": mean(depth_at_mixed),
            # `expectation - law(mean)` in metres of regolith, per cell, before
            # the catena factors. One-signed by the curvature of a single
            # saturating function: mixing first can only make the profile
            # thinner. A gap that read zero everywhere would say the correction
            # is inert on this map, which is a finding and not a reason to drop
            # the operator.
            "jensen_gap_m": {
                "land_mean": mean(regolith_gap),
                "min": float(regolith_gap[land].min()),
                "max": float(regolith_gap[land].max()),
                "percentiles": {str(p): float(np.percentile(regolith_gap[land], p))
                                for p in (50, 75, 95, 99)},
                "land_area_fraction_above_bar": float(np.average(
                    (np.abs(regolith_gap) > pedo["regolith"]["minimum_depth_m"]
                     )[land].astype(float), weights=lw)),
            },
            # The catena factors are cell-scale statistics and scale both orders
            # by the same number, so the clip that follows them is the only step
            # after the law where the two orders can part again.
            "post_catena_clipped_land_cells": int(np.count_nonzero(land
                                                                   & catena_clipped)),
            "transfer": regolith_ledger,
        },
        "land_means": {
            "weathering_intensity": mean(intensity),
            "sand": mean(texture["sand"]),
            "silt": mean(texture["silt"]),
            "clay": mean(texture["clay"]),
            "quartz_inert": mean(texture["quartz"]),
            "ph": mean(ph),
            "regolith_depth_m": mean(depth),
            "organic_fraction": mean(organic_fraction),
            "bulk_density_kg_m3": mean(bulk_density),
            "water_capacity_mm": mean(water_capacity),
            "water_capacity_mm_bracket": [mean(water_capacity_low),
                                          mean(water_capacity_high)],
            "column_base_m": column_base,
            "land_cells_regolith_below_column_base": int(
                np.count_nonzero(land & (depth > column_base))),
            "land_cells": int(np.count_nonzero(land)),
            "median_share_of_regolith_below_column_base": float(np.median(
                (1.0 - column_base / depth)[land & (depth > column_base)]))
                if np.any(land & (depth > column_base)) else 0.0,
            "bedrock_water_fraction": mean(bedrock_fraction),
            "soil_carbon_kg_m2": mean(carbon),
            "runoff_mm_per_earth_year": mean(runoff),
            "precipitation_mm_per_earth_year": mean(precip),
            "temperature_c": mean(temperature),
            "endorheic_fraction": mean(endorheic),
            "tan_slope": mean(tan_beta),
            "frost_cycling_fraction_of_orbit": mean(frost_fraction),
            "slope_fines_lost": mean(fines_loss),
        },
        "column_base_note": (
            "the plant-available capacity is cut at the land column property "
            "contract's declared base. Nothing in this pipeline draws water "
            "from below it: the physical column is fifteen layers to that "
            "depth, LPJ-GUESS's root distribution spans exactly those layers, "
            "and the contract's rootable base is shallower still. Water in "
            "regolith deeper than the base belongs to the aquifer term, which "
            "hydrography/config/land_water_ledger.yaml registers as the "
            "unrepresentable `transient_saturated_storage` under PLHY-4. "
            "WORLD-7702."),
        "regolith_note": (
            "LPJ-GUESS has a fixed 1.5 m physical profile. VesperInput consumes "
            "regolith depth and weathered-bedrock fraction by scaling each "
            "layer's texture-derived water capacity, but it does not change "
            "layer or root geometry; see pedology/README.md and the soil/land-"
            "surface hydraulic consistency audit."),
        "output": rel(output),
        "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=PROJECT_ROOT).stdout.strip() or None,
    }
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    report_path = ANALYSIS / "soil_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")

    means = report["land_means"]
    print(f"land cells          {len(rows)}")
    bracket = report["weathering_bracket"]
    print(f"weathering W        {means['weathering_intensity']:.3f} "
          f"(1.0 = Earth land mean), driver {selector}")
    print(f"  bracket           {bracket['by_runoff']:.3f} by runoff, "
          f"{bracket['by_precipitation']:.3f} by precipitation")
    print(f"  runoff ratio      {bracket['runoff_ratio']:.3f} "
          f"against Earth land's ~0.35")
    print(f"texture             sand {means['sand']:.3f}  silt {means['silt']:.3f}  "
          f"clay {means['clay']:.3f}   inert quartz {means['quartz_inert']:.3f}")
    print(f"pH                  {means['ph']:.2f}")
    print(f"regolith depth      {means['regolith_depth_m']:.2f} m "
          f"(catena: slope {means['tan_slope']:.3f}, frost cycling "
          f"{means['frost_cycling_fraction_of_orbit']:.3f} of the orbit)")
    print(f"organic fraction    {means['organic_fraction']:.4f} "
          f"from {means['soil_carbon_kg_m2']:.2f} kgC/m2")
    print(f"bulk density        {means['bulk_density_kg_m3']:.0f} kg/m3")
    lo, hi = means["water_capacity_mm_bracket"]
    print(f"water capacity      {means['water_capacity_mm']:.1f} mm, "
          f"{lo:.1f} to {hi:.1f} across the declared bracket "
          f"(ExoPlaSim's uniform default is 500)")
    print(f"bedrock water       {means['bedrock_water_fraction']:.3f} of soil "
          f"capacity per unit volume")
    print(f"\nwrote {report['output']}")
    print(f"      {rel(report_path)}")


if __name__ == "__main__":
    main()
