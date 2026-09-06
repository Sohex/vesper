#!/usr/bin/env python3
"""Compare matched Vesper and stock-4.1.1 soil-N transformation runs."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

from _paths import CONFIG, COMPONENT_ROOT, RUNS
from gridding import label_row_weights
from rungs import model_grid

ANALYSIS = COMPONENT_ROOT / "analysis"
WINDOW_YEARS = 10
QUANTITIES = (
    ("gross_nitrification", "soil_nflux.out", "GROSS_NITRIF"),
    ("net_nitrification", "soil_nflux.out", "NET_NITRIF"),
    ("gross_denitrification", "soil_nflux.out", "GROSS_DENITRIF"),
    ("net_denitrification", "soil_nflux.out", "NET_DENITRIF"),
    ("soil_nh3", "soil_nflux.out", "NH3"),
    ("soil_no", "soil_nflux.out", "NO"),
    ("soil_n2o", "soil_nflux.out", "N2O"),
    ("soil_n2", "soil_nflux.out", "N2"),
    ("plant_mineral_n_uptake", "nuptake.out", "Total"),
)


def load_manifest(run: Path) -> dict:
    path = run / "run_manifest.json"
    if not path.is_file():
        raise SystemExit(f"{path} is missing")
    return json.loads(path.read_text(encoding="utf-8"))


def table(path: Path, field: str) -> dict[tuple[float, float, int], float]:
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split()
    if field not in header:
        raise SystemExit(f"{path} has no {field} column")
    index = header.index(field)
    rows = {}
    for line in lines[1:]:
        values = line.split()
        key = (round(float(values[0]), 5), round(float(values[1]), 5),
               int(float(values[2])))
        rows[key] = float(values[index])
    return rows


def require_matched(vesper: dict, stock: dict) -> list[str]:
    checks = []

    def same(name: str, left, right) -> None:
        if left != right:
            raise SystemExit(f"arms are not matched on {name}: {left!r} != {right!r}")
        checks.append(name)

    if vesper["ntransform_comparison"]["profile"] != "vesper":
        raise SystemExit("first run is not the Vesper ntransform profile")
    if stock["ntransform_comparison"]["profile"] != "stock-4.1.1":
        raise SystemExit("second run is not the stock-4.1.1 ntransform profile")
    stock_binary = stock["ntransform_comparison"].get("binary_provenance") or {}
    if stock_binary.get("contract_version") != "vesper-ntransform-comparison-binary/1":
        raise SystemExit("stock binary lacks comparison-arm provenance")

    for key in ("nyear", "npatch", "root_seed", "ranks", "nfix_a", "nfix_b",
                "ifbvoc", "run_peatland", "ifmethane"):
        same(f"physical.{key}", vesper["physical"].get(key), stock["physical"].get(key))
    for key in ("driver", "soilmap", "pfts", "vesper_h"):
        same(f"inputs.{key}.sha256", vesper["inputs"][key]["sha256"],
             stock["inputs"][key]["sha256"])
    for key in ("forcing", "stochastic_randomness", "config_sha256", "source_build",
                "soil_driver_climate_contract"):
        same(key, vesper.get(key), stock.get(key))
    same("cells_simulated", vesper["cells_simulated"], stock["cells_simulated"])
    return checks


def summarize(vesper_rows: dict, stock_rows: dict) -> dict:
    if set(vesper_rows) != set(stock_rows):
        raise SystemExit("the two tables do not contain identical cell-year keys")
    years = sorted({key[2] for key in vesper_rows})
    if len(years) < WINDOW_YEARS:
        raise SystemExit(f"comparison needs {WINDOW_YEARS} retained years")
    selected = set(years[-WINDOW_YEARS:])
    cells = sorted({key[:2] for key in vesper_rows})
    v_cell = []
    s_cell = []
    for lon, lat in cells:
        v_cell.append(np.mean([vesper_rows[(lon, lat, year)] for year in selected]))
        s_cell.append(np.mean([stock_rows[(lon, lat, year)] for year in selected]))
    v = np.asarray(v_cell)
    s = np.asarray(s_cell)
    # The table is keyed by the latitude each cell carries rather than by an
    # index, so the weight is looked up on the model's rows and REFUSES a
    # latitude that is not one of them. `lib/gridding.py` owns both halves.
    _, nlat, _ = model_grid(yaml.safe_load(CONFIG.read_text(encoding="utf-8")))
    w = label_row_weights([lat for _, lat in cells], nlat,
                          what=f"the {len(cells)} simulated cells")
    delta = v - s
    vmean = float(np.average(v, weights=w))
    smean = float(np.average(s, weights=w))
    scale = max(abs(smean), 1.0e-30)
    return {
        "unit": "kgN/ha/simulation-year",
        "stock_area_weighted_mean": smean,
        "vesper_area_weighted_mean": vmean,
        "vesper_minus_stock": vmean - smean,
        "relative_change_percent": 100.0 * (vmean - smean) / scale,
        "cell_delta_quantiles": {
            "p05": float(np.quantile(delta, 0.05)),
            "p50": float(np.quantile(delta, 0.50)),
            "p95": float(np.quantile(delta, 0.95)),
        },
        "cell_sign_fraction": {
            "lower": float(np.mean(delta < 0.0)),
            "equal": float(np.mean(delta == 0.0)),
            "higher": float(np.mean(delta > 0.0)),
        },
        "cells": len(cells),
        "window_years": years[-WINDOW_YEARS:],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vesper_run")
    parser.add_argument("stock_run")
    args = parser.parse_args()
    vesper_run = RUNS / args.vesper_run
    stock_run = RUNS / args.stock_run
    vesper_manifest = load_manifest(vesper_run)
    stock_manifest = load_manifest(stock_run)
    invariants = require_matched(vesper_manifest, stock_manifest)

    results = {}
    cache = {}
    for name, filename, field in QUANTITIES:
        for arm, run in (("vesper", vesper_run), ("stock", stock_run)):
            key = (arm, filename, field)
            cache[key] = table(run / filename, field)
        results[name] = summarize(cache[("vesper", filename, field)],
                                  cache[("stock", filename, field)])

    report = {
        "contract_version": "vesper-ntransform-comparison/1",
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "vesper_run": args.vesper_run,
        "stock_run": args.stock_run,
        "comparison": "Vesper minus stock LPJ-GUESS 4.1.1",
        "window_complete_forcing_cycles": WINDOW_YEARS,
        "matched_invariants": invariants,
        "profiles": {
            "vesper": vesper_manifest["ntransform_comparison"],
            "stock": stock_manifest["ntransform_comparison"],
        },
        "quantities": results,
    }
    target = ANALYSIS / f"ntransform_{args.vesper_run}_{args.stock_run}.json"
    target.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(target)
    for name, values in results.items():
        print(f"{name:28s} stock={values['stock_area_weighted_mean']:.8g} "
              f"vesper={values['vesper_area_weighted_mean']:.8g} "
              f"change={values['relative_change_percent']:+.3f}%")


if __name__ == "__main__":
    main()
