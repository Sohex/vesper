"""Score an LPJ-GUESS run against the pre-registered productivity prediction.

`notes/productivity-prediction.md` fixed ten numbered predictions, each with a
hit band and a clear-miss band, before the model had run a single Vesper
gridcell. This computes what the run actually did and reports hit or miss per
line. It does not adjust anything: a miss is the result.

    python biosphere/scripts/score_prediction.py biosphere/runs/<run_id>

The bands below are transcribed from that document and are not to be edited to
fit a result. If a band turns out to be badly posed, say so in the note and
re-register; do not quietly widen it here.
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

from _paths import CONFIG, PROJECT_ROOT, climatology_path

import orbit

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = COMPONENT_ROOT / "analysis"

# Earth reference, pinned in the prediction so the comparison cannot drift.
EARTH_LAND_AREA_KM2 = 149.0e6
EARTH_TOTAL_NPP_PGC = (55.0, 60.0)
EARTH_LAND_MEAN_NPP = (370.0, 400.0)

GRASS_PFTS = ("C3G", "C4G")
BOREAL_PFTS = ("BNE", "BINE", "BNS")

# Cover above which a PFT counts as present in a cell. The registration said
# "tree FPC > 0.1" for prediction 5 and left 9 and 10 looser; the same threshold
# is applied to all three and stated here rather than chosen per line.
PRESENCE_FPC = 0.1


def read_table(path: Path) -> tuple[list[str], dict[tuple[float, float], list[float]]]:
    """Last simulated year per cell, from an LPJ-GUESS .out file."""
    lines = path.read_text().splitlines()
    header = lines[0].split()
    names = header[3:]
    latest: dict[tuple[float, float], tuple[int, list[float]]] = {}
    for line in lines[1:]:
        parts = line.split()
        if len(parts) < 3 + len(names):
            continue
        key = (round(float(parts[0]), 2), round(float(parts[1]), 2))
        year = int(float(parts[2]))
        if key not in latest or year > latest[key][0]:
            latest[key] = (year, [float(v) for v in parts[3:3 + len(names)]])
    return names, {k: v for k, (_, v) in latest.items()}


def band(value: float, hit: tuple[float, float],
         miss: tuple[float, float]) -> str:
    if hit[0] <= value <= hit[1]:
        return "HIT"
    if miss[0] <= value <= miss[1]:
        return "outside hit band, inside clear-miss band"
    return "CLEAR MISS"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, help="an LPJ-GUESS run directory")
    parser.add_argument("--climatology", type=Path, default=None,
                        help="the climatology the run was driven from, for the "
                             "forcing-derived lines and cell areas")
    args = parser.parse_args()

    run = args.run
    manifest_path = run / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {}
    config = yaml.safe_load(CONFIG.read_text())
    climatology = args.climatology or climatology_path()

    # Simulation years are this world's years. Everything registered is per Earth
    # year. Getting this wrong looks like missing every line low by 2.
    to_earth = manifest.get("simulation_years_to_earth_years")
    if to_earth is None:
        to_earth = 365.2568983 / orbit.orbital_year_days(config)

    with nc.Dataset(climatology) as data:
        lat = np.asarray(data["lat"][:], dtype=float)
        lon = np.asarray(data["lon"][:], dtype=float)
        land = np.asarray(data["lsm"][0], dtype=float) > 0.5
        temperature = np.asarray(data["tas"][:], dtype=float).mean(axis=0) - 273.15
        precip = (np.asarray(data["pr"][:], dtype=float).mean(axis=0)
                  * 1000.0 * 86400.0 * 365.2425)

    radius_km = 6371.0 * float(config["planet"]["radius_earth"])
    weight = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, len(lon)))
    cell_km2 = weight / weight.sum() * 4.0 * np.pi * radius_km ** 2
    lon_signed = np.where(lon > 180.0, lon - 360.0, lon)

    npp_names, npp = read_table(run / "anpp.out")
    lai_names, lai = read_table(run / "lai.out")
    fpc_names, fpc = read_table(run / "fpc.out")

    # Map every simulated cell onto the grid so areas can weight it.
    simulated = np.zeros(land.shape, dtype=bool)
    npp_grid = np.zeros(land.shape)
    lai_grid = np.zeros(land.shape)
    tree_grid = np.zeros(land.shape)
    grass_grid = np.zeros(land.shape)
    c4_grid = np.zeros(land.shape)
    boreal_grid = np.zeros(land.shape)

    tree_idx = [i for i, n in enumerate(fpc_names)
                if n not in GRASS_PFTS and n != "Total"]
    grass_idx = [i for i, n in enumerate(fpc_names) if n in GRASS_PFTS]
    c4_idx = [i for i, n in enumerate(fpc_names) if n == "C4G"]
    boreal_idx = [i for i, n in enumerate(fpc_names) if n in BOREAL_PFTS]
    npp_total = npp_names.index("Total")
    lai_total = lai_names.index("Total")

    for j in range(len(lat)):
        for i in range(len(lon)):
            key = (round(float(lon_signed[i]), 2), round(float(lat[j]), 2))
            row = fpc.get(key)
            if row is None:
                continue
            simulated[j, i] = True
            npp_grid[j, i] = npp[key][npp_total]
            lai_grid[j, i] = lai[key][lai_total]
            tree_grid[j, i] = sum(row[k] for k in tree_idx)
            grass_grid[j, i] = sum(row[k] for k in grass_idx)
            c4_grid[j, i] = sum(row[k] for k in c4_idx)
            boreal_grid[j, i] = sum(row[k] for k in boreal_idx)

    covered = simulated & land
    n_land = int(land.sum())
    n_done = int(covered.sum())
    if n_done == 0:
        raise SystemExit("no simulated cells matched the climatology grid")
    partial = n_done < n_land

    area = cell_km2[covered]
    land_area_km2 = float(cell_km2[land].sum())

    def wmean(field: np.ndarray) -> float:
        return float(np.average(field[covered], weights=area))

    def wfrac(mask: np.ndarray) -> float:
        return float(np.average(mask[covered].astype(float), weights=area))

    # gC/m2 per Earth year, from kgC/m2 per simulation year.
    land_mean_npp = wmean(npp_grid) * 1000.0 * to_earth
    total_npp_pgc = land_mean_npp * land_area_km2 * 1e6 / 1e15
    earth_land_mean = sum(EARTH_LAND_MEAN_NPP) / 2.0
    earth_total = sum(EARTH_TOTAL_NPP_PGC) / 2.0

    # Prediction 7 is a property of the forcing, not of LPJ-GUESS: the
    # registration defined it through the Miami model, so it is recomputed the
    # same way rather than invented from model output.
    def miami(t, p):
        return np.minimum(3000 / (1 + np.exp(1.315 - 0.119 * t)),
                          3000 * (1 - np.exp(-0.000664 * p))) * 0.45
    temp_limit = 3000 / (1 + np.exp(1.315 - 0.119 * temperature))
    prec_limit = 3000 * (1 - np.exp(-0.000664 * precip))
    water_limited = wfrac(prec_limit < temp_limit)

    results = [
        ("1  land-mean NPP relative to Earth",
         land_mean_npp / earth_land_mean, (0.7, 1.4), (0.5, 2.0), "x"),
        ("2  land-mean NPP, gC/m2 per Earth year",
         land_mean_npp, (280.0, 560.0), (200.0, 800.0), ""),
        ("3  total NPP relative to Earth",
         total_npp_pgc / earth_total, (1.6, 2.8), (1.2, 3.5), "x"),
        ("4  total NPP, PgC per Earth year",
         total_npp_pgc, (90.0, 175.0), (60.0, 240.0), ""),
        ("5  tree cover, land fraction with tree FPC > 0.1",
         wfrac(tree_grid > PRESENCE_FPC), (0.40, 0.70), (0.25, 0.80), ""),
        ("6  effectively barren land, LAI < 0.5",
         wfrac(lai_grid < 0.5), (0.08, 0.25), (0.0, 0.40), ""),
        ("7  water-limited land (from the forcing, not the run)",
         water_limited, (0.50, 0.70), (0.35, 1.0), ""),
        ("8  grass cover exceeds tree cover",
         1.0 if wmean(grass_grid) > wmean(tree_grid) else 0.0,
         (1.0, 1.0), (1.0, 1.0), " (1 = yes)"),
        ("9  C4 grasses above 10% of land",
         wfrac(c4_grid > PRESENCE_FPC), (0.10, 1.0), (0.05, 1.0), ""),
        ("10 boreal needleleaf below 10% of land",
         wfrac(boreal_grid > PRESENCE_FPC), (0.0, 0.10), (0.0, 0.20), ""),
    ]

    print(f"run              {run.name}")
    print(f"cells scored     {n_done} of {n_land} land cells"
          + ("   PARTIAL RUN, treat every line as indicative" if partial else ""))
    print(f"simulation year  {to_earth:.4f} Earth years")
    print(f"land area        {land_area_km2 / 1e6:.1f}e6 km2\n")
    scored = []
    for name, value, hit, miss, unit in results:
        verdict = band(value, hit, miss)
        print(f"  {name:52s} {value:8.3f}{unit:9s} {verdict}")
        scored.append({"prediction": name, "value": value, "hit_band": list(hit),
                       "clear_miss_band": list(miss), "verdict": verdict})

    hits = sum(1 for r in scored if r["verdict"] == "HIT")
    print(f"\n  {hits} of {len(scored)} predictions hit")

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "run": str(run),
        "run_manifest": manifest,
        "climatology": str(climatology.relative_to(PROJECT_ROOT)),
        "cells_scored": n_done,
        "land_cells": n_land,
        "partial": partial,
        "simulation_years_to_earth_years": to_earth,
        "land_area_km2": land_area_km2,
        "land_mean_npp_gc_m2_earth_year": land_mean_npp,
        "total_npp_pgc_earth_year": total_npp_pgc,
        "earth_reference": {
            "land_area_km2": EARTH_LAND_AREA_KM2,
            "total_npp_pgc": list(EARTH_TOTAL_NPP_PGC),
            "land_mean_npp": list(EARTH_LAND_MEAN_NPP),
        },
        "predictions": scored,
        "hits": hits,
        "note": ("Bands are transcribed from notes/productivity-prediction.md and "
                 "are not to be edited to fit a result. A miss is the result."),
        "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=PROJECT_ROOT).stdout.strip() or None,
    }
    out_dir = ANALYSIS / run.name
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "prediction_score.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwrote {path.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
