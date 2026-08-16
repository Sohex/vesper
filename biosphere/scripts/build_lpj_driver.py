"""Build the binary driver file LPJ-GUESS reads, from climate and lithology.

One self-describing file carrying the gridlist, a soil code per cell derived from
World Orogen lithology, and binned climate. `vesperinput.cpp` reads it.

Three things are worth knowing about what goes in.

**Land comes from `surface_class`, never `land_mask`.** The two disagree by 1.9%
of the planet, all of it dry closed-basin floor below sea level that `land_mask`
would flood. That is the whole point of the Orogen fork.

**Insolation is supplied as net downward surface shortwave.** `driver.cpp`'s
`NETSWRAD_TS` path then applies no albedo correction of its own, which is what we
want: this project computes surface albedo from lithology and ExoPlaSim has
already used it, so the number is better than driver.cpp's global 0.17 constant.

**Bin order defines the calendar.** Bin 0 becomes day 0 of the simulation year,
and `build_vesper_header.py` fits the declination phase on that assumption. The
two must be regenerated together after any orbit change.

**One climatology or many.** Pass several to `--climatology` and each becomes one
year of forcing, in the order given; `vesperinput` cycles through them, so the
spin-up sees the whole sequence rather than one arbitrary phase of it. One file
is a fixed climate. Several are how a stellar cycle reaches the biosphere, since
a single repeating year cannot represent a variable star at all.

    python biosphere/scripts/build_lpj_driver.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import struct
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import netCDF4 as nc
import numpy as np
import yaml

from _paths import CONFIG, GENERATED, PROJECT_ROOT, climatology_path

import builds
import orbit
from gridding import land_fraction_of_class
from orogen import Export

MAGIC = b"VESPDRV4"   # V2 regolith depth, V3 bedrock water, V4 multiple years

# Coordinate precision shared with pedology/scripts/build_soil.py, so the soil
# map keys match exactly. See where lon_signed is rounded.
COORD_DECIMALS = 4
PROVENANCE_BYTES = 64
KELVIN = 273.15

# LPJ soil codes, from the table in modules/soilinput.h:
#   0 ice, 1 coarse, 2 medium, 3 fine, 4 medium-coarse, 5 fine-coarse,
#   6 fine-medium, 7 fine-medium-coarse, 8 organic, 9 vertisols
#
# Mapping Orogen's rock classes onto them is a judgement call, so it is written
# down rather than buried. The reasoning is weathering product, not parent rock
# hardness: what matters to a plant is the texture of the regolith the rock
# breaks down into and the water it can hold.
#
# Nothing maps to 8 (organic), because organic soils are a product of the
# biosphere we are about to model and cannot be an input to it. Nothing maps to
# 0 (ice) either; glaciated ground is handled by the climate model, not here.
SOIL_CODE_BY_ROCK = {
    "water":               2,  # unused, land only, but keep the table total
    "morb":                4,  # basalt weathers to clay-rich but stony ground
    "oib":                 4,
    "flood_basalt":        6,  # deep basalt saprolite, finer than arc material
    "arc_basalt":          4,
    "arc_andesite":        4,
    "rift_bimodal":        4,
    "granite":             1,  # granite grus is coarse and sandy
    "granodiorite":        1,
    "gneiss":              4,  # gneiss weathers coarse but with more clay
    "schist":              6,  # phyllosilicates give fine, platy regolith
    "quartzite":           1,  # almost pure quartz sand
    "melange":             7,  # tectonic mixture, so a mixed texture
    "shelf_clastic":       7,  # interbedded sandstone and shale
    "carbonate":           3,  # karst residuum is clay-rich terra rossa
    "foreland_clastic":    7,
    "continental_clastic": 5,
    "pelagic":             3,  # abyssal clay
    "evaporite":           9,  # vertisol: the shrink-swell salt-affected case
    "playa_clastic":       3,  # playa mud, fine and poorly drained
}


def project_relative(path: Path) -> str:
    """Relative to the project when it is inside it, absolute when it is not.

    Provenance should read cleanly for tracked inputs without crashing on a
    scratch path outside the tree, which is exactly what /tmp climatologies are
    during a test.
    """
    path = Path(path).resolve()
    return str(path.relative_to(PROJECT_ROOT)
               if path.is_relative_to(PROJECT_ROOT) else path)


def area_weights(lat: np.ndarray, nlon: int) -> np.ndarray:
    w = np.cos(np.deg2rad(lat))
    return (w / w.sum())[:, None] * np.ones((1, nlon)) / nlon


def soil_codes(config: dict, land: np.ndarray) -> tuple[np.ndarray, dict]:
    """Dominant soil code per grid cell, integrated from the native mesh.

    Not sampled from the gridded export. That export resamples categorical fields
    by the region containing the cell centre, which at T42 throws away most of
    what is in a cell. Instead each soil code's share of a cell's *land* area is
    accumulated over the 2.5M-region mesh and the largest share wins.

    Dominant rather than averaged, because a soil code is categorical: the mean
    of code 1 and code 3 is code 2, which is not what half granite and half
    carbonate behaves like.
    """
    mesh = Export(builds.mesh_export(config))
    grid_dir = builds.grid_export(config, str(config["model"]["resolution"]).upper())
    rock = mesh.surface_rock
    classes = {c["id"]: c["code"] for c in mesh.manifest["lithology"]["rockClasses"]}

    unknown = sorted({classes[i] for i in np.unique(rock) if i in classes}
                     - set(SOIL_CODE_BY_ROCK))
    if unknown:
        raise SystemExit(
            f"rock classes with no soil-code mapping: {unknown}. Add them to "
            f"SOIL_CODE_BY_ROCK and record why."
        )

    # Accumulate land-area share per soil code, then take the argmax per cell.
    shares: dict[int, np.ndarray] = {}
    for rock_id, code in classes.items():
        soil = SOIL_CODE_BY_ROCK[code]
        selected = rock == rock_id
        if not selected.any():
            continue
        share = land_fraction_of_class(mesh, grid_dir, selected)
        shares[soil] = shares.get(soil, np.zeros_like(share)) + share

    codes_present = sorted(shares)
    stacked = np.stack([shares[c] for c in codes_present])
    dominant = np.array(codes_present, dtype=np.int32)[np.argmax(stacked, axis=0)]

    # Cells with no land in the mesh but flagged land by the climatology's mask
    # would otherwise take whichever code sorts first. Medium is the neutral
    # default and the count is reported rather than hidden.
    empty = stacked.sum(axis=0) <= 0.0
    dominant[empty] = 2
    stranded = int(np.sum(empty & land))

    total = max(int(land.sum()), 1)
    summary = {
        "method": "dominant land-area share, integrated over the native mesh",
        "cells_with_no_mesh_land_defaulted_to_medium": stranded,
        "soil_code_share_of_land_cells": {
            int(c): round(float(np.sum(dominant[land] == c) / total), 4)
            for c in sorted(np.unique(dominant[land]))
        },
    }
    return dominant, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--climatology", type=Path, nargs="+", default=None,
                        help="one or more climatologies, each one year of "
                             "forcing, cycled in the order given")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--soil-map", type=Path, default=None,
                        help="pedology soilmap.txt, for its regolith depth "
                             "column. Without it every cell is given the full "
                             "profile depth, which is LPJ-GUESS's own default.")
    parser.add_argument("--ndep", type=float, default=0.5,
                        help="nitrogen deposition, kgN/ha/yr. A declared "
                             "assumption: this world has no deposition field and "
                             "no industry. Default is a low pre-industrial-like "
                             "value; report the sensitivity, do not tune it.")
    args = parser.parse_args()

    config = yaml.safe_load(CONFIG.read_text())
    climatologies = list(args.climatology) if args.climatology else [climatology_path()]
    for path in climatologies:
        if not Path(path).is_file():
            raise SystemExit(f"{path} does not exist")
    climatology = Path(climatologies[0])

    year_length = int(round(orbit.orbital_year_days(config)))
    co2_ppm = float(config["atmosphere"]["pCO2_bar"]) / 1.0 * 1e6

    # Each climatology contributes one year, stacked as [year][bin][lat][lon].
    tas_y, pr_y, rss_y, dtr_y = [], [], [], []
    lat = lon = lsm = None
    for path in climatologies:
        with nc.Dataset(path) as data:
            this_lat = np.asarray(data["lat"][:], dtype=float)
            this_lon = np.asarray(data["lon"][:], dtype=float)
            if lat is None:
                lat, lon = this_lat, this_lon
                lsm = np.asarray(data["lsm"][0], dtype=float)
            elif not (np.allclose(lat, this_lat) and np.allclose(lon, this_lon)):
                raise SystemExit(
                    f"{path} is on a different grid from {climatologies[0]}; "
                    f"every year has to share one grid")
            tas_y.append(np.asarray(data["tas"][:], dtype=float) - KELVIN)
            pr_y.append(np.asarray(data["pr"][:], dtype=float) * 1000.0 * 86400.0)
            rss_y.append(np.asarray(data["rss"][:], dtype=float))
            dtr_y.append(np.maximum(
                np.asarray(data["maxt"][:], dtype=float)
                - np.asarray(data["mint"][:], dtype=float), 0.0))
    tas = np.stack(tas_y)   # [year][bin][lat][lon]
    pr = np.stack(pr_y)
    rss = np.stack(rss_y)
    dtr = np.stack(dtr_y)
    nyears = tas.shape[0]

    nbins = tas.shape[1]
    if nbins != 12:
        raise SystemExit(f"driver format expects 12 bins per year, got {nbins}")

    land = lsm > 0.5
    codes, soil_summary = soil_codes(config, land)

    # Bin length in days, from the same month lengths the patched Date uses, so
    # a precipitation total is the total for exactly the days it is spread over.
    base = year_length // 12
    bin_days = np.array([base] * 12, dtype=float)
    bin_days[-1] += year_length - base * 12

    # ExoPlaSim's longitudes run 0..360; LPJ-GUESS expects -180..180.
    #
    # Rounded to COORD_DECIMALS because LPJ-GUESS's soil map lookup keys on
    # std::pair<double,double> and compares it exactly. pedology/build_soil.py
    # writes its coordinates to the same precision, so the two round-trip to
    # identical doubles and the lookup needs no search radius. That agreement is
    # a contract between the two scripts, not a coincidence.
    lon_signed = np.round(np.where(lon > 180.0, lon - 360.0, lon), COORD_DECIMALS)
    lat = np.round(lat, COORD_DECIMALS)

    provenance = f"{config.get('source_build')}|{climatology.parent.name}".encode()
    provenance = provenance[:PROVENANCE_BYTES].ljust(PROVENANCE_BYTES, b"\0")

    output = args.output or (GENERATED / "vesper_driver.bin")
    output.parent.mkdir(parents=True, exist_ok=True)

    # Regolith depth, keyed on the same rounded coordinates the soil map uses.
    # Absent, every cell gets LPJ-GUESS's full 1.5 m profile, which is what the
    # unpatched model assumes anyway.
    default_depth_m = 1.5
    depth_by_coord: dict[tuple[float, float], float] = {}
    bedrock_by_coord: dict[tuple[float, float], float] = {}
    # Without a soil map, sub-bedrock layers keep the fresh-rock minimum. It is a
    # declared physical value, not a numerical guard; see pedogenesis.yaml.
    default_bedrock_fraction = 0.05
    soil_map = args.soil_map
    if soil_map is None:
        candidate = PROJECT_ROOT / "pedology" / "data" / "soilmap.txt"
        soil_map = candidate if candidate.is_file() else None
    if soil_map is not None:
        lines = soil_map.read_text().splitlines()
        header = lines[0].split()
        for needed in ("depth", "bedrockfrac"):
            if needed not in header:
                raise SystemExit(
                    f"{soil_map} has no {needed} column; rebuild it with "
                    f"build_soil.py")
        depth_column = header.index("depth")
        bedrock_column = header.index("bedrockfrac")
        for line in lines[1:]:
            parts = line.split()
            key = (round(float(parts[0]), COORD_DECIMALS),
                   round(float(parts[1]), COORD_DECIMALS))
            depth_by_coord[key] = float(parts[depth_column])
            bedrock_by_coord[key] = float(parts[bedrock_column])
        print(f"regolith depth from {soil_map.name}: {len(depth_by_coord)} cells")
    else:
        print("no soil map found; every cell gets the full 1.5 m profile")

    rows = np.argwhere(land)
    missing_depth = 0
    with output.open("wb") as handle:
        handle.write(MAGIC)
        handle.write(struct.pack("<iiii", len(rows), nbins, year_length, nyears))
        handle.write(struct.pack("<dd", co2_ppm, args.ndep))
        handle.write(provenance)
        for j, i in rows:
            handle.write(struct.pack("<dd", float(lon_signed[i]), float(lat[j])))
            handle.write(struct.pack("<ii", int(codes[j, i]), 0))
            key = (float(lon_signed[i]), float(lat[j]))
            depth = depth_by_coord.get(key)
            if depth is None:
                depth = default_depth_m
                missing_depth += 1
            handle.write(struct.pack("<d", depth))
            handle.write(struct.pack(
                "<d", bedrock_by_coord.get(key, default_bedrock_fraction)))
            # Flattened [year][bin], matching what vesperinput indexes.
            handle.write(tas[:, :, j, i].astype("<f8").tobytes())
            handle.write((pr[:, :, j, i] * bin_days[None, :]).astype("<f8").tobytes())
            handle.write(rss[:, :, j, i].astype("<f8").tobytes())
            handle.write(dtr[:, :, j, i].astype("<f8").tobytes())

    weights = area_weights(lat, len(lon))
    lw = weights[land]
    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "climatologies": [project_relative(p) for p in climatologies],
        "climatology_sha256": [hashlib.sha256(Path(p).read_bytes()).hexdigest()
                               for p in climatologies],
        "years_of_climate": nyears,
        "years_note": ("Cycled by vesperinput, so spin-up sees the whole "
                       "sequence. One year is a fixed climate; several are how a "
                       "stellar cycle reaches the biosphere."),
        "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
        "source_build": config.get("source_build"),
        "stellar_spectrum": config.get("radiation", {}).get("stellar_spectrum"),
        "flux_earth": float(config["orbit"]["baseline_flux_earth"]),
        "orbital_year_earth_days": orbit.orbital_year_days(config),
        "year_length_days": year_length,
        "bins_per_year": nbins,
        "bin_days": bin_days.tolist(),
        "land_cells": int(len(rows)),
        "regolith_depth_source": (project_relative(soil_map)
                                  if soil_map else "none, full profile assumed"),
        "cells_without_depth": missing_depth,
        "co2_ppm": co2_ppm,
        "ndep_kgn_ha_yr": args.ndep,
        "ndep_note": (
            "Declared, not measured. This world has no deposition field and no "
            "industry, so any value is an assumption and NPP inherits it."),
        "insolation": "NETSWRAD_TS, net downward surface shortwave (rss), W/m2",
        "land_definition": "lsm from the climatology, itself built from surface_class",
        "land_mean_temperature_c": float(
            np.average(tas.mean(axis=(0, 1))[land], weights=lw)),
        "land_mean_precip_mm_per_earth_year": float(
            np.average(pr.mean(axis=(0, 1))[land], weights=lw) * 365.2425),
        "land_mean_net_sw_w_m2": float(
            np.average(rss.mean(axis=(0, 1))[land], weights=lw)),
        "soil": soil_summary,
        "soil_code_mapping": SOIL_CODE_BY_ROCK,
        "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=PROJECT_ROOT).stdout.strip() or None,
    }
    report_path = output.with_name(output.stem + "_provenance.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n")

    print(f"land cells     {len(rows)}")
    print(f"year length    {year_length} days, bins {bin_days.astype(int).tolist()}")
    print(f"climate years  {nyears} "
          f"({'cycled' if nyears > 1 else 'fixed climate, repeated'})")
    print(f"CO2            {co2_ppm:.0f} ppm     N deposition {args.ndep} kgN/ha/yr")
    print(f"land means     {report['land_mean_temperature_c']:.2f} C, "
          f"{report['land_mean_precip_mm_per_earth_year']:.0f} mm/Earth-yr, "
          f"{report['land_mean_net_sw_w_m2']:.1f} W/m2 net SW")
    print("soil codes     " + ", ".join(
        f"{k}:{v:.0%}" for k, v in soil_summary["soil_code_share_of_land_cells"].items()))
    print(f"\nwrote {project_relative(output)} "
          f"({output.stat().st_size / 1e6:.1f} MB)")
    print(f"      {report_path.name}")


if __name__ == "__main__":
    main()
