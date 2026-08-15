#!/usr/bin/env python3
"""Average a contiguous post-equilibrium ExoPlaSim segment into climatologies."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re

from netCDF4 import Dataset
import numpy as np

from _paths import ANALYSIS


COORDINATES = {"time", "lat", "lon", "lev", "levp", "fourier", "modes"}
ORBITAL_GEOMETRY = {"nu", "lambda", "zdec", "rdist", "rasc"}


def average_files(paths: list[Path], output: Path, product: str) -> None:
    """Write a compressed, same-grid mean across corresponding model orbits."""
    if not paths:
        raise ValueError("No input files supplied")
    for path in paths:
        if not path.is_file():
            raise RuntimeError(f"Missing climatology input {path}")

    output.parent.mkdir(parents=True, exist_ok=True)
    with Dataset(paths[0]) as template, Dataset(output, "w", format="NETCDF4") as dst:
        for name, dim in template.dimensions.items():
            dst.createDimension(name, len(dim))

        dst.setncatts({name: template.getncattr(name) for name in template.ncattrs()})
        dst.setncattr("climatology_product", product)
        dst.setncattr("climatology_orbit_count", len(paths))
        dst.setncattr("climatology_start_year_index", int(paths[0].name.split(".")[-2]))
        dst.setncattr("climatology_end_year_index", int(paths[-1].name.split(".")[-2]))
        dst.setncattr("climatology_created_utc", datetime.now(timezone.utc).isoformat())
        dst.setncattr("climatology_source_files", ",".join(path.name for path in paths))

        for name, source in template.variables.items():
            fill = source.getncattr("_FillValue") if "_FillValue" in source.ncattrs() else None
            kwargs = {"zlib": True, "complevel": 4, "shuffle": True}
            if fill is not None:
                kwargs["fill_value"] = fill
            target = dst.createVariable(name, source.dtype, source.dimensions, **kwargs)
            target.setncatts(
                {
                    attr: source.getncattr(attr)
                    for attr in source.ncattrs()
                    if attr != "_FillValue"
                }
            )

            first = np.ma.asarray(source[:])
            should_average = "time" in source.dimensions and name not in (
                COORDINATES | ORBITAL_GEOMETRY
            )
            if should_average:
                total = np.ma.asarray(first, dtype=np.float64)
                for path in paths[1:]:
                    with Dataset(path) as nc:
                        if name not in nc.variables or nc[name].shape != source.shape:
                            raise RuntimeError(f"Variable mismatch for {name} in {path}")
                        total += np.ma.asarray(nc[name][:], dtype=np.float64)
                target[:] = total / len(paths)
            else:
                for path in paths[1:]:
                    with Dataset(path) as nc:
                        if name not in nc.variables or nc[name].shape != source.shape:
                            raise RuntimeError(f"Variable mismatch for {name} in {path}")
                        if not np.ma.allclose(first, nc[name][:], rtol=1e-6, atol=1e-7):
                            raise RuntimeError(
                                f"Coordinate/orbital variable {name} differs in {path}"
                            )
                target[:] = first


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--start-year", type=int, required=True)
    parser.add_argument("--end-year", type=int, required=True)
    parser.add_argument("--output", type=Path, default=ANALYSIS / "climatology")
    parser.add_argument("--label", default="baseline")
    args = parser.parse_args()
    if args.end_year < args.start_year:
        raise ValueError("--end-year must not precede --start-year")
    if not re.fullmatch(r"[a-z0-9_-]+", args.label):
        raise ValueError("--label may contain only lowercase letters, digits, _ and -")

    run_dir = args.run_dir.resolve()
    years = range(args.start_year, args.end_year + 1)
    regular = [run_dir / f"MOST.{year:05d}.nc" for year in years]
    snapshots = [run_dir / "snapshots" / f"MOST_SNAP.{year:05d}.nc" for year in years]
    output_dir = args.output.resolve()
    regular_output = output_dir / f"{args.label}_regular_climatology.nc"
    snapshot_output = output_dir / f"{args.label}_snapshot_climatology.nc"
    average_files(regular, regular_output, "time-bin means averaged across model orbits")
    average_files(snapshots, snapshot_output, "instantaneous orbital snapshots averaged across model orbits")

    manifest_path = run_dir / "run_manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest.setdefault("climatologies", {})[args.label] = {
            "start_year_index": args.start_year,
            "end_year_index": args.end_year,
            "orbit_count": args.end_year - args.start_year + 1,
            "regular": str(regular_output),
            "snapshots": str(snapshot_output),
        }
        # Retain the original convenience key for existing tooling.
        if args.label == "baseline":
            manifest["climatology"] = manifest["climatologies"][args.label]
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(
        json.dumps(
            {"regular": str(regular_output), "snapshots": str(snapshot_output)},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
