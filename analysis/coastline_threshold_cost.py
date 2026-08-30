#!/usr/bin/env python3
"""Measure SPAT-5's area/store coastline cost on every exported ladder rung."""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

SCRIPT = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT.parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "lib"))
sys.path.insert(0, str(PROJECT_ROOT / "exoplasim/scripts"))

import builds
from orogen import Export, INLAND_WATER
from build_boundary_conditions import build

CONFIG = PROJECT_ROOT / "config/planet.yaml"
OUTPUT = PROJECT_ROOT / "analysis/coastline_threshold_cost.json"
LADDER = ("T21", "T42", "T85", "T127", "T170")


def main() -> None:
    config = yaml.safe_load(CONFIG.read_text())
    mesh = Export(builds.mesh_export(config))
    threshold = float(config["model"]["geography_land_threshold"])
    gravity = float(config["planet"]["gravity_m_s2"])
    rows = []
    for rung in LADDER:
        grid = builds.grid_export(config, rung)
        if not grid.is_dir():
            raise SystemExit(f"{grid} is absent; SPAT-5 requires the whole ladder")
        result = build(mesh, grid, threshold, gravity)
        ledger = result["coastline_ledger"]
        rows.append({"rung": rung, "grid": str(grid.relative_to(PROJECT_ROOT)),
                     "cells_below_mesh_resolution": result[
                         "cells_below_mesh_resolution"],
                     "ledger": ledger})
        print(f"{rung}: dropped land "
              f"{100 * ledger['land_dropped_to_ocean']['of_mesh_land']:.3f}%, "
              f"promoted water "
              f"{100 * ledger['water_promoted_to_land']['of_model_land']:.3f}%, "
              f"volume {100 * ledger['land_volume_closure_relative']:+.3f}%")

    inland_area_fraction = float(
        mesh.cell_area[mesh.surface_class == INLAND_WATER].sum()
        / mesh.cell_area.sum())
    report = {
        "contract_version": "vesper-partial-surface-area-store/1",
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source_build": config["source_build"],
        "terrain_hash": mesh.terrain_hash,
        "threshold": threshold,
        "rungs": rows,
        "third_class": {
            "mesh_area_fraction": inland_area_fraction,
            "status": ("absent_not_small" if inland_area_fraction == 0.0
                       else "present"),
            "detail": ("the active export contains exactly zero inherited "
                       "INLAND_WATER area; solved lakes are a later land "
                       "partition and must not be folded into this zero"),
        },
        "selection": {
            "status": "awaiting_accepted_baseline_flux_bracket",
            "tile_model_selected": False,
            "model_form_bracket_selected": False,
            "reason": ("author decision: cost bracket first, choose a tile "
                       "model only if the measured response is material"),
        },
        "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
    }
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n")
    print(f"wrote {OUTPUT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
