#!/usr/bin/env python3
"""Validate and record a completed segment after post-run bookkeeping failed."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import yaml

from _paths import CONFIG

from continue_exoplasim import validate_year, year_diagnostics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--start-year", type=int, required=True)
    parser.add_argument("--end-year", type=int, required=True)
    parser.add_argument("--seasonal-output", action="store_true")
    parser.add_argument("--config", type=Path, default=CONFIG)
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    model = config["model"]
    diagnostics = []
    for year in range(args.start_year, args.end_year + 1):
        validate_year(
            run_dir / f"MOST.{year:05d}.nc",
            expected_times=int(model["regular_output_bins_per_orbit"]),
        )
        if args.seasonal_output:
            validate_year(
                run_dir / "snapshots" / f"MOST_SNAP.{year:05d}.nc",
                expected_times=int(model["seasonal_samples_per_orbit"]),
            )
        diagnostics.append(year_diagnostics(run_dir / f"MOST.{year:05d}.nc"))

    manifest_path = run_dir / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    existing = manifest.setdefault("segments", [])
    if any(
        segment["start_year_index"] == args.start_year
        and segment["end_year_index"] == args.end_year
        for segment in existing
    ):
        raise RuntimeError("Segment is already recorded in the manifest")
    is_climatology = bool(
        args.seasonal_output
        and "equilibrium_cutoff_year_index" in manifest
        and args.start_year > int(manifest["equilibrium_cutoff_year_index"])
    )
    existing.append(
        {
            "start_year_index": args.start_year,
            "end_year_index": args.end_year,
            "seasonal_output": args.seasonal_output,
            "purpose": "post_equilibrium_climatology" if is_climatology else "spinup",
            "recovered_after_bookkeeping_failure_utc": datetime.now(timezone.utc).isoformat(),
            "diagnostics": diagnostics,
        }
    )
    manifest["completed_orbits"] = args.end_year + 1
    manifest["status"] = "climatology_complete" if is_climatology else "spinup_in_progress"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()
