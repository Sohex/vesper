#!/usr/bin/env python3
"""Generate world_state.json: what this project currently knows about Vesper.

Generated, not hand-maintained. Every prose summary in this project has drifted
from its artifacts at least once, usually within a commit of the number changing,
so the state file reads the artifacts instead of restating them. Anything here
that cannot be derived is in `curated`, which lives in this file under version
control rather than being edited into the JSON.

    python scripts/world_state.py

Re-run it after anything that changes a build, a run, or a verdict.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
SCHEMA_VERSION = 1

# --- things no artifact can tell us -----------------------------------------
# Judgement, not data. Kept in the generator so it is reviewed in diffs.
CURATED = {
    "open_questions": [
        "Does the climate response follow the albedo arithmetic? Land-mean albedo "
        "on carved-zoned is 0.1722 against the baseline's 0.2240, predicting about "
        "299.4 K. The chain is four linear steps across a response we have "
        "repeatedly declined to extrapolate.",
        "The habitable flux depends on the biosphere, which is assumed rather than "
        "modelled. Bare rock and vegetated reach 290-293 K at non-overlapping "
        "fluxes, so LPJ-GUESS determines not just the climate but where the planet "
        "should be placed.",
        "Iteration 1 systematically over-carved: carving removes bright evaporite, "
        "the world warms, and marginal basins would have stayed closed. The "
        "overshoot should be measured on iteration 2 by re-evaluating the carved "
        "set against the new climate.",
        "Roughly 228 basins with zero catchment runoff were preserved but have "
        "P > E over open water, so a lake there grows and spills. They should go "
        "through the marginal band with a retain fraction, not be flipped to 0.",
    ],
    "known_uncertainties": [
        {"item": "playa_clastic albedo", "value": 0.30, "range": [0.25, 0.33],
         "worth_w_m2": 1.99, "worth_k": 1.87,
         "note": "Dominant remaining albedo lever: 18.8% of land against the salt "
                 "crust's 2.1%, 2.48 W/m2 per 0.10. Itself a mixture of clay playa "
                 "and desert-varnished fan gravel. Recorded, not tuned."},
        {"item": "salt crust extent", "value": 0.106, "range": [0.10, 0.40],
         "note": "Degenerate with crust albedo in the effective value: 10.6% at "
                 "0.50 and 20% at 0.40 are the same planet from the climate side. "
                 "Only an independent measure of flooding frequency separates "
                 "them, and the lake solver is not validated enough to provide it."},
        {"item": "lake levels", "value": None,
         "note": "The equilibrium solver has never been validated against anything, "
                 "unlike Penman which was checked against the model's own ocean "
                 "cells to 1.7%. Least-checked numbers in the pipeline."},
        {"item": "energy budget closure", "value": 0.45,
         "note": "Net TOA runs 0.27-0.45 W/m2 more negative than net surface flux "
                 "across every run, which cannot both be true in steady state. The "
                 "|mean TOA| < 0.5 criterion is being applied inside a bias of its "
                 "own size, so margins of 0.004 either way are not resolvable."},
        {"item": "retain fraction", "value": None,
         "note": "Measures distance from threshold, not incision capacity. A basin "
                 "that barely trickles over its sill is cut the same as one that "
                 "pours. Should become a function of overflow discharge."},
    ],
    "next_steps": [
        "Finish the carved-zoned baseline at 0.96 and compare against the 299.4 K "
        "prediction.",
        "Choose the flux from what it measures, not from the estimate.",
        "Supply equilibrium lake levels as --water-level so basin floors stop being "
        "exported as dry ground.",
        "Run LPJ-GUESS on the resulting climatology and replace the assumed "
        "biosphere.",
        "Iteration-2 carve verdict, with retain as a function of discharge and the "
        "228 zero-runoff spillers handled.",
        "One T85 equilibrium once terrain and biosphere settle, for regional "
        "products.",
        "Stellar cycle last, 0.91-1.01 S-Earth over 8 Earth years.",
    ],
}


def git_commit() -> str | None:
    try:
        return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return None


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else None


def builds() -> dict:
    out = {}
    for d in sorted((ROOT / "source").iterdir()):
        man = d / "exoplasim-T42" / "manifest.json"
        if not man.is_file():
            continue
        m = json.loads(man.read_text(encoding="utf-8"))
        b, lit = m["basins"], m["lithology"]
        comp = {x["code"]: x["fraction"] for x in lit["compositionLand"]}
        out[d.name] = {
            "terrain_hash": m["hashes"]["finalElevation"],
            "catalogue_hash": m["hashes"]["basinCatalogue"],
            "seed": m["seed"],
            "regions": m["numRegions"],
            "selection_source": b.get("selectionSource"),
            "drainage_hypothesis": {k: v for k, v in
                                    (b.get("drainageHypothesis") or {}).items()
                                    if k != "entries"},
            "preserved_basins": b["counts"]["preserved"],
            "land_fraction": m.get("landSeaMask", {}).get("landFractionBySurfaceClass"),
            "endorheic_fraction_of_land":
                (b.get("drainageConsistency") or {}).get("fractionOfLand"),
            "lithology_land_fractions": {k: round(v, 5) for k, v in
                                         sorted(comp.items(), key=lambda kv: -kv[1])[:6]},
            "grids": sorted(p.name.replace("exoplasim-", "") for p in d.iterdir()
                            if p.is_dir() and p.name.startswith("exoplasim-")),
            "superseded": (d / "SUPERSEDED.md").is_file(),
        }
    return out


def climate_runs() -> list:
    rows = []
    for man in sorted((ROOT / "exoplasim" / "runs").glob("*/run_manifest.json")):
        d = json.loads(man.read_text(encoding="utf-8"))
        run = man.parent
        conv = d.get("convergence_assessment") or {}
        # The stellar-cycle runs use fixed_orbit_parameters rather than
        # derived_parameters, and the earliest runs predate several keys, so read
        # everything defensively: a state file that crashes on an old artifact is
        # worse than one that reports a null.
        derived = d.get("derived_parameters") or d.get("fixed_orbit_parameters") or {}
        cfg = d.get("source_config") or {}
        model = cfg.get("model") or {}
        rows.append({
            "run_id": d.get("run_id", run.name),
            "status": d.get("status"),
            "orbits_complete": len(sorted(run.glob("MOST.*.nc"))),
            "flux_earth": derived.get("stellar_flux_ratio_earth"),
            "cycle": d.get("cycle", {}).get("case") if isinstance(d.get("cycle"), dict) else None,
            "resolution": model.get("resolution"),
            "land_albedo_source": model.get("land_albedo_source"),
            "source_build": cfg.get("source_build"),
            "glaciers": (cfg.get("surface") or {}).get("glaciers", {}).get("enabled"),
            "converged": conv.get("pass"),
            "mean_surface_temperature_k": (conv.get("metrics") or {}).get(
                "temperature_mean_k") or _mean_ts(run),
        })
    return rows


def _mean_ts(run: Path, last: int = 5) -> float | None:
    """Global mean surface temperature over the final orbits, from the output.

    Computed here rather than taken from a convergence assessment, because a run
    has output long before anything assesses it, and a state file reporting null
    for a finished run is worse than one that does the arithmetic.
    """
    files = sorted(run.glob("MOST.*.nc"))
    if not files:
        return None
    try:
        from netCDF4 import Dataset
        from numpy.polynomial.legendre import leggauss
        vals = []
        for f in files[-last:]:
            with Dataset(f) as ds:
                w = leggauss(ds.dimensions["lat"].size)[1][::-1]
                a = np.asarray(ds["ts"][:])
                while a.ndim > 2:
                    a = a.mean(axis=0)
                vals.append(float((a.mean(axis=1) * w).sum() / w.sum()))
        return round(float(np.mean(vals)), 3)
    except Exception:
        return None


def hydrography() -> dict:
    rep = read_json(ROOT / "hydrography" / "data" / "hydrography_report.json") or {}
    verdict = read_json(ROOT / "hydrography" / "analysis" / "carve_verdict.json") or {}
    carve = read_json(ROOT / "hydrography" / "data" / "carve_list.json") or {}
    return {
        "terrain_hash": (rep.get("source") or {}).get("terrain_hash"),
        "drainage": rep.get("drainage"),
        "basins": rep.get("basins"),
        "carve_verdict": {
            "bounds": verdict.get("bounds"),
            "agreement": verdict.get("agreement"),
            "penman_ocean_validation": verdict.get("penman_ocean_validation"),
        },
        "carve_list": {
            "counts": carve.get("counts"),
            "terrain_hash": carve.get("terrain_hash"),
            "climate": carve.get("climate"),
        } if carve else None,
    }


def current_climate() -> dict | None:
    reports = sorted((ROOT / "exoplasim" / "analysis").glob("*/baseline_climate_report.json"))
    if not reports:
        return None
    latest = reports[-1]
    d = json.loads(latest.read_text(encoding="utf-8"))
    koppen = d.get("koppen_land_area_fractions") or {}
    biome = d.get("biome_land_area_fractions") or {}
    return {
        "source": str(latest.relative_to(ROOT)),
        "global_metrics": d.get("global_metrics"),
        "koppen_top": dict(sorted(koppen.items(), key=lambda kv: -kv[1])[:6]),
        "biome_top": dict(sorted(biome.items(), key=lambda kv: -kv[1])[:6]),
        "caveat": "Computed on precarve-unzoned with evaporite as one class at "
                  "0.50. Both have since been superseded. Quote as the 0.96 "
                  "pre-carve baseline, never as Vesper's climate.",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=ROOT / "world_state.json")
    args = ap.parse_args()

    config = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    state = {
        "schema_version": SCHEMA_VERSION,
        "generated": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit(),
        "generator": "scripts/world_state.py",
        "note": "Generated from the artifacts. Do not edit; re-run the generator. "
                "Judgement that no artifact can supply is in `curated`, which "
                "lives in the generator under version control.",
        "planet": config["planet"],
        "star": config["star"],
        "orbit": config["orbit"],
        "atmosphere": config["atmosphere"],
        "active_configuration": {
            "source_build": config.get("source_build"),
            "model": config["model"],
            "surface": config["surface"],
            "radiation": config["radiation"],
            "stellar_cycle": config.get("stellar_cycle"),
        },
        "builds": builds(),
        "climate_runs": climate_runs(),
        "current_climate": current_climate(),
        "hydrography": hydrography(),
        "curated": CURATED,
    }
    args.output.write_text(json.dumps(state, indent=2, default=str) + "\n",
                           encoding="utf-8")

    b = state["builds"]
    print(f"builds        : {len(b)}  (active: {config.get('source_build')})")
    print(f"climate runs  : {len(state['climate_runs'])}, "
          f"{sum(1 for r in state['climate_runs'] if r['converged'])} converged")
    hy = state["hydrography"]
    if hy.get("basins"):
        print(f"basins        : {hy['basins']['count']}")
    if hy.get("carve_list") and hy["carve_list"].get("counts"):
        print(f"carve list    : {hy['carve_list']['counts']}")
    print(f"open questions: {len(CURATED['open_questions'])}, "
          f"uncertainties: {len(CURATED['known_uncertainties'])}")
    print(f"\nwrote {args.output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
