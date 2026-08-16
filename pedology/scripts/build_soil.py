"""Weather the lithology into soil, under this world's climate.

Between the climate model and the vegetation model there is a step neither of
them does: turning rock into something a plant can root in. Parent material sets
what minerals are available; climate sets how far they have been converted;
relief and erosion set how much regolith survives; and the biosphere sets the
organic fraction, which is what makes this a loop rather than a stage.

Emits the soil map LPJ-GUESS's `SoilInput` reads, one row per land cell:

    Lon Lat sand clay silt orgc ph bulkdensity cn soilc

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

from _paths import (ANALYSIS, CONFIG, DATA, PEDOGENESIS, PROJECT_ROOT,
                    climatology_path)

import builds
import orbit
from gridding import land_fraction_of_class
from orogen import Export

# Coordinate precision shared with biosphere/scripts/build_lpj_driver.py.
# LPJ-GUESS keys its soil map on an exactly-compared pair of doubles, so both
# files must round identically or every lookup misses.
COORD_DECIMALS = 4

EARTH_YEAR_DAYS = 365.2425
KELVIN = 273.15

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
    rock = mesh.surface_rock
    classes = {c["id"]: c["code"] for c in mesh.manifest["lithology"]["rockClasses"]}

    fractions: dict[str, np.ndarray] = {}
    for rock_id, code in classes.items():
        selected = rock == rock_id
        if selected.any():
            fractions[code] = land_fraction_of_class(mesh, grid_dir, selected)
    return fractions, mesh, grid_dir


def weathering_intensity(runoff_mm_yr: np.ndarray, temperature_c: np.ndarray,
                         params: dict) -> np.ndarray:
    """Walker-Hays-Kasting weathering intensity, normalised to Earth land means.

    Dimensionless and clipped. It is an intensity rather than a rate: the time
    integral is folded into the reference, because soil age is not known on
    either world well enough to carry explicitly.
    """
    q_ratio = np.maximum(runoff_mm_yr, 0.0) / params["reference_runoff_mm_per_earth_year"]
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

    weatherable = np.clip(1.0 - quartz - clay, 0.0, 1.0)
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


def regolith_depth(intensity: np.ndarray, relief_m: np.ndarray,
                   runoff_mm_yr: np.ndarray, erodibility: np.ndarray,
                   params: dict, weathering_ref: float) -> np.ndarray:
    """Steady-state regolith thickness from production against erosion.

    Heimsath's exponential soil production function balanced against an erosion
    rate proportional to erodibility, local relief and runoff. Where erosion
    outruns production the depth floors: that is bare rock.
    """
    erosion = (erodibility
               * np.maximum(relief_m, 0.0) / params["erosion_reference_relief_m"]
               * np.maximum(runoff_mm_yr, 0.0) / weathering_ref)
    production = intensity
    with np.errstate(divide="ignore", invalid="ignore"):
        depth = params["e_folding_depth_m"] * np.log(
            np.where(erosion > 1e-9, production / erosion, np.inf))
    depth = np.where(np.isfinite(depth), depth, params["maximum_depth_m"])
    return np.clip(depth, params["minimum_depth_m"], params["maximum_depth_m"])


def soil_ph(fractions: dict[str, np.ndarray], runoff_mm_yr: np.ndarray,
            endorheic: np.ndarray, params: dict, reference_runoff: float
            ) -> np.ndarray:
    """Parent pH, leached down by drainage, pushed up where drainage is closed."""
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
    climatology = args.climatology or climatology_path()
    if not climatology.is_file():
        raise SystemExit(f"{climatology} does not exist")

    with nc.Dataset(climatology) as data:
        lat = np.asarray(data["lat"][:], dtype=float)
        lon = np.asarray(data["lon"][:], dtype=float)
        land = np.asarray(data["lsm"][0], dtype=float) > 0.5
        temperature = np.asarray(data["tas"][:], dtype=float).mean(axis=0) - KELVIN
        # Runoff and precipitation are rates; annualise on Earth years so the
        # Earth-calibrated weathering law is fed the units it was fitted in.
        runoff = np.asarray(data["mrro"][:], dtype=float).mean(axis=0) \
            * 1000.0 * 86400.0 * EARTH_YEAR_DAYS
        precip = np.asarray(data["pr"][:], dtype=float).mean(axis=0) \
            * 1000.0 * 86400.0 * EARTH_YEAR_DAYS
        elevation = np.asarray(data["sg"][:], dtype=float).mean(axis=0) \
            / float(config["planet"]["gravity_m_s2"])

    fractions, mesh, grid_dir = lithology_fractions(config)

    # Local relief, as the spread of surface height across each cell's neighbours.
    # A stand-in for slope that needs no extra field and no mesh gradient.
    relief = np.zeros_like(elevation)
    for axis in (0, 1):
        forward = np.roll(elevation, 1, axis=axis)
        backward = np.roll(elevation, -1, axis=axis)
        relief = np.maximum(relief, np.maximum(np.abs(elevation - forward),
                                               np.abs(elevation - backward)))

    erodibility = np.zeros_like(elevation)
    total_share = np.zeros_like(elevation)
    for entry in mesh.manifest["lithology"]["rockClasses"]:
        share = fractions.get(entry["code"])
        if share is None:
            continue
        erodibility += share * float(entry["erodibility"])
        total_share += share
    covered = total_share > 1e-9
    erodibility[covered] /= total_share[covered]
    erodibility[~covered] = 1.0

    endorheic = land_fraction_of_class(mesh, grid_dir, mesh.is_endorheic.astype(bool))

    # Both moisture drivers, so the bracket is visible whichever one is selected.
    intensity_by_runoff = weathering_intensity(runoff, temperature, pedo["weathering"])
    intensity_by_precip = weathering_intensity(precip, temperature, pedo["weathering"])
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
    depth = regolith_depth(intensity, relief, runoff, erodibility,
                           pedo["regolith"],
                           pedo["weathering"]["reference_runoff_mm_per_earth_year"])
    ph = soil_ph(fractions, runoff, endorheic, pedo["ph"],
                 pedo["weathering"]["reference_runoff_mm_per_earth_year"])

    if args.soil_carbon:
        carbon = read_soil_carbon(args.soil_carbon, lon, lat)
    else:
        carbon = np.full(temperature.shape,
                         float(pedo["organic"]["initial_soil_carbon_kg_m2"]))
    organic_fraction, bulk_density = organic_properties(carbon, pedo["organic"])

    # Plant-available water capacity, mm: volumetric capacity from texture times
    # the depth of regolith that actually exists. This is what ExoPlaSim's dwmax
    # bucket should be, and its runoff is literally the overflow of that bucket,
    # so it is the field that closes the loop back to the climate.
    water = pedo["water"]
    volumetric = (water["volumetric_capacity_by_texture"]["sand"] * texture["sand"]
                  + water["volumetric_capacity_by_texture"]["silt"] * texture["silt"]
                  + water["volumetric_capacity_by_texture"]["clay"] * texture["clay"]
                  + water["volumetric_capacity_organic"] * organic_fraction)
    water_capacity = np.clip(volumetric * depth * 1000.0,
                             water["minimum_mm"], water["maximum_mm"])

    # Plant-available water below the bedrock contact, as a fraction of what the
    # soil above holds per unit volume. A function of weathering, because that
    # is what turns impermeable rock into saprock and then saprolite. See the
    # anchors in pedogenesis.yaml; this is not a small correction.
    bedrock = pedo["regolith"]["bedrock_water"]
    bedrock_fraction = (bedrock["minimum"]
                        + (bedrock["maximum"] - bedrock["minimum"])
                        * (1.0 - np.exp(-bedrock["shape"] * intensity)))

    DATA.mkdir(parents=True, exist_ok=True)
    # Resolved so a relative or out-of-tree --output does not break the
    # provenance record's relative_to(PROJECT_ROOT).
    output = (args.output or (DATA / "soilmap.txt")).resolve()
    lon_signed = np.where(lon > 180.0, lon - 360.0, lon)
    rows = np.argwhere(land)
    with output.open("w") as handle:
        # `depth` and `awc` are extra columns LPJ-GUESS's own SoilInput ignores
        # (it matches columns by name and skips what it does not know). They are
        # here so one artifact carries the whole soil: biosphere reads depth to
        # scale water capacity, and exoplasim reads awc to set its bucket.
        handle.write("Lon Lat sand clay silt orgc ph bulkdensity cn soilc "
                     "depth awc bedrockfrac\n")
        for j, i in rows:
            handle.write(
                f"{lon_signed[i]:.{COORD_DECIMALS}f} {lat[j]:.{COORD_DECIMALS}f} "
                f"{texture['sand'][j, i]:.4f} {texture['clay'][j, i]:.4f} "
                f"{texture['silt'][j, i]:.4f} {organic_fraction[j, i]:.5f} "
                f"{ph[j, i]:.3f} {bulk_density[j, i]:.1f} "
                f"{pedo['organic']['carbon_nitrogen_ratio']:.1f} "
                f"{carbon[j, i]:.4f} "
                f"{depth[j, i]:.4f} {water_capacity[j, i]:.2f} "
                f"{bedrock_fraction[j, i]:.4f}\n")

    weights = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, len(lon)))
    lw = weights[land]

    def mean(field: np.ndarray) -> float:
        return float(np.average(field[land], weights=lw))

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "iteration": args.iteration,
        "closes_loop_with": (str(args.soil_carbon) if args.soil_carbon
                             else "nothing; iteration 0 has no biosphere"),
        "climatology": str(climatology.relative_to(PROJECT_ROOT)),
        "climatology_sha256": hashlib.sha256(climatology.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
        "pedogenesis_sha256": hashlib.sha256(PEDOGENESIS.read_bytes()).hexdigest(),
        "source_build": config.get("source_build"),
        "terrain_hash": mesh.terrain_hash,
        "orbital_year_earth_days": orbit.orbital_year_days(config),
        "land_cells": int(len(rows)),
        "moisture_variable": selector,
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
            "bedrock_water_fraction": mean(bedrock_fraction),
            "soil_carbon_kg_m2": mean(carbon),
            "runoff_mm_per_earth_year": mean(runoff),
            "precipitation_mm_per_earth_year": mean(precip),
            "temperature_c": mean(temperature),
            "endorheic_fraction": mean(endorheic),
        },
        "regolith_note": (
            "LPJ-GUESS 4.1.1 has a fixed 1.5 m profile and does not consume "
            "regolith depth. It is computed and reported so the gap is visible; "
            "see notes/model.md."),
        "output": (str(output.relative_to(PROJECT_ROOT))
                   if output.is_relative_to(PROJECT_ROOT) else str(output)),
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
    print(f"regolith depth      {means['regolith_depth_m']:.2f} m (not consumed)")
    print(f"organic fraction    {means['organic_fraction']:.4f} "
          f"from {means['soil_carbon_kg_m2']:.2f} kgC/m2")
    print(f"bulk density        {means['bulk_density_kg_m3']:.0f} kg/m3")
    print(f"water capacity      {means['water_capacity_mm']:.1f} mm "
          f"(ExoPlaSim's uniform default is 500)")
    print(f"bedrock water       {means['bedrock_water_fraction']:.3f} of soil "
          f"capacity per unit volume")
    print(f"\nwrote {report['output']}")
    print(f"      {report_path.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
