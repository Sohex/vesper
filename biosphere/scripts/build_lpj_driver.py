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

MAGIC = b"VESPDRV1"
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
    parser.add_argument("--climatology", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--ndep", type=float, default=0.5,
                        help="nitrogen deposition, kgN/ha/yr. A declared "
                             "assumption: this world has no deposition field and "
                             "no industry. Default is a low pre-industrial-like "
                             "value; report the sensitivity, do not tune it.")
    args = parser.parse_args()

    config = yaml.safe_load(CONFIG.read_text())
    climatology = args.climatology or climatology_path()
    if not climatology.is_file():
        raise SystemExit(f"{climatology} does not exist")

    year_length = int(round(orbit.orbital_year_days(config)))
    co2_ppm = float(config["atmosphere"]["pCO2_bar"]) / 1.0 * 1e6

    with nc.Dataset(climatology) as data:
        lat = np.asarray(data["lat"][:], dtype=float)
        lon = np.asarray(data["lon"][:], dtype=float)
        lsm = np.asarray(data["lsm"][0], dtype=float)
        tas = np.asarray(data["tas"][:], dtype=float) - KELVIN
        pr = np.asarray(data["pr"][:], dtype=float) * 1000.0 * 86400.0  # mm/day
        rss = np.asarray(data["rss"][:], dtype=float)                   # W/m2
        maxt = np.asarray(data["maxt"][:], dtype=float)
        mint = np.asarray(data["mint"][:], dtype=float)

    nbins = tas.shape[0]
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
    lon_signed = np.where(lon > 180.0, lon - 360.0, lon)

    provenance = f"{config.get('source_build')}|{climatology.parent.name}".encode()
    provenance = provenance[:PROVENANCE_BYTES].ljust(PROVENANCE_BYTES, b"\0")

    output = args.output or (GENERATED / "vesper_driver.bin")
    output.parent.mkdir(parents=True, exist_ok=True)

    rows = np.argwhere(land)
    with output.open("wb") as handle:
        handle.write(MAGIC)
        handle.write(struct.pack("<iiii", len(rows), nbins, year_length, 0))
        handle.write(struct.pack("<dd", co2_ppm, args.ndep))
        handle.write(provenance)
        for j, i in rows:
            handle.write(struct.pack("<dd", float(lon_signed[i]), float(lat[j])))
            handle.write(struct.pack("<ii", int(codes[j, i]), 0))
            handle.write(tas[:, j, i].astype("<f8").tobytes())
            handle.write((pr[:, j, i] * bin_days).astype("<f8").tobytes())
            handle.write(rss[:, j, i].astype("<f8").tobytes())
            handle.write(np.maximum(maxt[:, j, i] - mint[:, j, i], 0.0)
                         .astype("<f8").tobytes())

    weights = area_weights(lat, len(lon))
    lw = weights[land]
    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "climatology": str(climatology.relative_to(PROJECT_ROOT)),
        "climatology_sha256": hashlib.sha256(climatology.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
        "source_build": config.get("source_build"),
        "stellar_spectrum": config.get("radiation", {}).get("stellar_spectrum"),
        "flux_earth": float(config["orbit"]["baseline_flux_earth"]),
        "orbital_year_earth_days": orbit.orbital_year_days(config),
        "year_length_days": year_length,
        "bins_per_year": nbins,
        "bin_days": bin_days.tolist(),
        "land_cells": int(len(rows)),
        "co2_ppm": co2_ppm,
        "ndep_kgn_ha_yr": args.ndep,
        "ndep_note": (
            "Declared, not measured. This world has no deposition field and no "
            "industry, so any value is an assumption and NPP inherits it."),
        "insolation": "NETSWRAD_TS, net downward surface shortwave (rss), W/m2",
        "land_definition": "lsm from the climatology, itself built from surface_class",
        "land_mean_temperature_c": float(np.average(tas.mean(axis=0)[land], weights=lw)),
        "land_mean_precip_mm_per_earth_year": float(
            np.average(pr.mean(axis=0)[land], weights=lw) * 365.2425),
        "land_mean_net_sw_w_m2": float(np.average(rss.mean(axis=0)[land], weights=lw)),
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
    print(f"CO2            {co2_ppm:.0f} ppm     N deposition {args.ndep} kgN/ha/yr")
    print(f"land means     {report['land_mean_temperature_c']:.2f} C, "
          f"{report['land_mean_precip_mm_per_earth_year']:.0f} mm/Earth-yr, "
          f"{report['land_mean_net_sw_w_m2']:.1f} W/m2 net SW")
    print("soil codes     " + ", ".join(
        f"{k}:{v:.0%}" for k, v in soil_summary["soil_code_share_of_land_cells"].items()))
    print(f"\nwrote {output.relative_to(PROJECT_ROOT)} "
          f"({output.stat().st_size / 1e6:.1f} MB)")
    print(f"      {report_path.name}")


if __name__ == "__main__":
    main()
