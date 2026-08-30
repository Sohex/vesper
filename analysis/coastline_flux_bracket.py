#!/usr/bin/env python3
"""Bracket SPAT-5's missing partial-cell flux on an accepted baseline.

This does not pretend that an offline mixture is a tiled atmosphere.  It holds
the accepted baseline climate fixed and replaces only the fraction that the
binary 0.5 mask assigns to the wrong surface class.  The unknown opposite-tile
rate is bounded by observed cells of that class in the same Gaussian row,
expanding symmetrically only when a row carries fewer than the declared minimum
number of references.  A result is consequently a model-form bracket, not a
prediction of the coupled response of separate tile temperatures and stores.

There is deliberately no fallback from the configured accepted baseline.  A
plausible number from an older terrain is worse than no number, which is the
state SPAT-5 is in until canonical-10m-carve2 has an accepted baseline.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml
from netCDF4 import Dataset

SCRIPT = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT.parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "lib"))
sys.path.insert(0, str(PROJECT_ROOT / "exoplasim/scripts"))

import builds
import climatology
import gridding
from build_boundary_conditions import build
from orogen import Export
from paths import climatology_path, rel

CONFIG = PROJECT_ROOT / "config/planet.yaml"
CONTRACT = PROJECT_ROOT / "config/partial_surface.yaml"
OUTPUT = PROJECT_ROOT / "analysis/coastline_flux_bracket.json"
RAW_FIELDS = ("evap", "hfss", "hfls", "rss", "rls")
MM_PER_DAY = 86400.0 * 1000.0


def reference_rows(class_mask: np.ndarray, row: int,
                   minimum: int) -> tuple[np.ndarray, int]:
    """Reference cells, expanding symmetrically by whole Gaussian rows."""
    nlat, _ = class_mask.shape
    for distance in range(nlat):
        lo, hi = max(0, row - distance), min(nlat, row + distance + 1)
        selected = np.zeros_like(class_mask, dtype=bool)
        selected[lo:hi] = class_mask[lo:hi]
        if int(selected.sum()) >= minimum:
            return selected, distance
    raise ValueError("a surface class has fewer references than the declared minimum")


def bracket_field(field: np.ndarray, land_fraction: np.ndarray,
                  binary_land: np.ndarray, cell_area: np.ndarray,
                  minimum_references: int) -> dict:
    """Global interval for fractional-minus-binary flux, with cell ledgers."""
    field = np.asarray(field, dtype=float)
    fraction = np.asarray(land_fraction, dtype=float)
    land = np.asarray(binary_land, dtype=bool)
    area = np.asarray(cell_area, dtype=float)
    if not (field.shape == fraction.shape == land.shape == area.shape):
        raise ValueError("field, land fraction, binary mask and area must share a shape")
    if not (np.isfinite(field).all() and np.isfinite(fraction).all()
            and np.isfinite(area).all()):
        raise ValueError("flux bracket inputs must be finite")
    if np.any((fraction < 0.0) | (fraction > 1.0)) or np.any(area <= 0.0):
        raise ValueError("land fractions must be in [0,1] and cell areas positive")

    partial = (fraction > 0.0) & (fraction < 1.0)
    lower = np.zeros_like(field)
    upper = np.zeros_like(field)
    expansions: list[dict] = []
    for row, col in np.argwhere(partial):
        target_land = not bool(land[row, col])
        target_class = land if target_land else ~land
        refs, distance = reference_rows(target_class, int(row), minimum_references)
        values = field[refs]
        missing_share = (fraction[row, col] if target_land
                         else 1.0 - fraction[row, col])
        candidates = missing_share * (values - field[row, col])
        lower[row, col] = float(candidates.min())
        upper[row, col] = float(candidates.max())
        expansions.append({
            "row": int(row), "column": int(col),
            "missing_surface": "land" if target_land else "ocean",
            "missing_share": float(missing_share),
            "row_expansion": int(distance),
            "reference_cells": int(refs.sum()),
        })

    lower_total = float(np.sum(lower * area, dtype=np.float64))
    upper_total = float(np.sum(upper * area, dtype=np.float64))
    sphere_area = float(np.sum(area, dtype=np.float64))
    return {
        "global_mean_delta_lower": lower_total / sphere_area,
        "global_mean_delta_upper": upper_total / sphere_area,
        "extensive_delta_lower": lower_total,
        "extensive_delta_upper": upper_total,
        "extensive_closure_residual": max(
            abs(lower_total - float(np.dot(lower.ravel(), area.ravel()))),
            abs(upper_total - float(np.dot(upper.ravel(), area.ravel())))),
        "partial_cells": int(partial.sum()),
        "reference_inventory": {
            "cells": expansions,
            "expanded_cells": int(sum(x["row_expansion"] > 0 for x in expansions)),
            "maximum_row_expansion": int(max(
                (x["row_expansion"] for x in expansions), default=0)),
            "minimum_references_observed": int(min(
                (x["reference_cells"] for x in expansions), default=0)),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    cfg = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    declaration = yaml.safe_load(CONTRACT.read_text(encoding="utf-8"))
    # The one accepted-baseline resolver. It raises on the current null instead
    # of selecting one of the stale files still present in the analysis tree.
    climate_path = climatology_path(root=PROJECT_ROOT)
    if not climate_path.is_file():
        raise SystemExit(f"accepted baseline does not exist: {rel(climate_path)}")

    mesh = Export(builds.mesh_export(cfg))
    grid_dir = builds.grid_export(cfg)
    built = build(mesh, grid_dir, float(cfg["model"]["geography_land_threshold"]),
                  float(cfg["planet"]["gravity_m_s2"]))
    fraction = np.asarray(built["land_fraction"], dtype=float)
    binary = np.asarray(built["land_mask"], dtype=bool)
    spec = gridding.gaussian_grid(*fraction.shape)
    radius = float(cfg["planet"]["radius_earth"]) * 6_371_000.0
    area = spec.cell_area(radius)
    minimum = int(declaration["flux_bracket"]["minimum_reference_cells"])

    with Dataset(climate_path) as ds:
        source_build = getattr(ds, "vesper_source_build", None)
        if source_build != cfg["source_build"]:
            raise SystemExit(
                f"{rel(climate_path)} belongs to source build {source_build!r}, "
                f"not active build {cfg['source_build']!r}")
        if int(getattr(ds, "low_io", 1)) != 0:
            raise SystemExit("accepted baseline must be built from NLOWIO=0 output")
        missing = [name for name in (*RAW_FIELDS, "lsm", "time")
                   if name not in ds.variables]
        if missing:
            raise SystemExit(f"accepted baseline lacks fields: {', '.join(missing)}")
        got_shape = (len(ds.dimensions["lat"]), len(ds.dimensions["lon"]))
        if got_shape != fraction.shape:
            raise SystemExit(f"baseline grid {got_shape} != export grid {fraction.shape}")
        lsm = climatology.annual_mean_of(ds, "lsm")
        if not np.array_equal(lsm >= 0.5, binary):
            disagree = int(np.count_nonzero((lsm >= 0.5) != binary))
            raise SystemExit(
                f"accepted baseline land mask differs from the active threshold "
                f"in {disagree} cells; array identity/orientation is not established")
        fields = {name: climatology.annual_mean_of(ds, name)
                  for name in RAW_FIELDS}
        fields["surface_energy"] = sum(fields[n] for n in ("rss", "rls", "hfss", "hfls"))
        run_id = getattr(ds, "vesper_run_id", None)

    results = {}
    for name, field in fields.items():
        row = bracket_field(field, fraction, binary, area, minimum)
        row["units"] = "m s-1" if name == "evap" else "W m-2"
        if name == "evap":
            row["global_mean_delta_lower_mm_day"] = (
                row["global_mean_delta_lower"] * MM_PER_DAY)
            row["global_mean_delta_upper_mm_day"] = (
                row["global_mean_delta_upper"] * MM_PER_DAY)
        results[name] = row

    report = {
        "contract_version": declaration["flux_bracket"]["contract_version"],
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": "measured_same_climate_model_form_bracket",
        "interpretation": declaration["flux_bracket"]["interpretation"],
        "source_build": cfg["source_build"],
        "terrain_hash": mesh.terrain_hash,
        "rung": cfg["model"]["resolution"],
        "threshold": float(cfg["model"]["geography_land_threshold"]),
        "accepted_baseline": rel(climate_path),
        "accepted_baseline_run_id": run_id,
        "accepted_baseline_sha256": hashlib.sha256(climate_path.read_bytes()).hexdigest(),
        "minimum_reference_cells": minimum,
        "method": declaration["flux_bracket"]["missing_tile_estimator"],
        "fields": results,
    }
    convergence_path = (PROJECT_ROOT / "exoplasim" / "analysis" / "convergence"
                        / f"{run_id}_convergence.json")
    if not convergence_path.is_file():
        raise SystemExit(
            f"accepted baseline has no convergence report: {rel(convergence_path)}")
    convergence = json.loads(convergence_path.read_text(encoding="utf-8"))
    storage_tolerance = float(
        convergence["criteria_provenance"]["storage_tolerance_w_m2"])
    surface_energy_bound = max(
        abs(results["surface_energy"]["global_mean_delta_lower"]),
        abs(results["surface_energy"]["global_mean_delta_upper"]),
    )
    tile_required = surface_energy_bound > storage_tolerance
    report["selection"] = {
        "status": ("tile_model_selected_implementation_pending" if tile_required
                   else "model_form_bracket_selected"),
        "materiality_test": declaration["selection"]["materiality_test"],
        "reference": rel(convergence_path),
        "state_storage_tolerance_w_m2": storage_tolerance,
        "maximum_absolute_surface_energy_bracket_w_m2": surface_energy_bound,
        "tile_model_selected": tile_required,
        "model_form_bracket_selected": not tile_required,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    for name, row in results.items():
        print(f"{name}: {row['global_mean_delta_lower']:+.6g} to "
              f"{row['global_mean_delta_upper']:+.6g} {row['units']}")
    print(f"wrote {rel(args.output)}")


if __name__ == "__main__":
    main()
