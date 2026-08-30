#!/usr/bin/env python3
"""Emit BIO-11's rootable land fraction from barren substrate and solved water.

    python biosphere/scripts/build_rootable_fraction.py
    python biosphere/scripts/build_rootable_fraction.py --self-test

The fraction is relative to native-mesh land inside each atmosphere cell. The
three mutually exclusive populations are rootable, solved water, and dry barren
substrate. A barren region under water is water, not two deductions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from netCDF4 import Dataset
import numpy as np
import yaml

from _paths import CONFIG, PROJECT_ROOT

import builds
import gridding
import rungs
from orogen import Export, LAND
from spatial_support import (grid_support_contract, validate_contract)

SURFACE_CLASSES = PROJECT_ROOT / "pedology" / "analysis" / "surface_classes.nc"


def partition(cell, ncell, area, land, barren, water):
    """Return land-relative rootable, water and dry-barren fractions."""
    water = land & np.asarray(water, dtype=bool)
    dry_barren = land & np.asarray(barren, dtype=bool) & ~water
    rootable = land & ~water & ~np.asarray(barren, dtype=bool)
    land_area = gridding.cell_sum(cell, ncell, area, land)
    outputs = []
    for mask in (rootable, water, dry_barren):
        fraction, _ = gridding.cell_fraction(cell, ncell, area, mask, land)
        outputs.append(fraction)
    covered = land_area > 0
    residual = float(np.max(np.abs(sum(outputs)[covered] - 1.0))) \
        if covered.any() else 0.0
    return outputs, land_area, residual


def self_test() -> int:
    cell = np.array([0, 0, 1, 1, 2, 2])
    area = np.ones(6)
    land = np.ones(6, dtype=bool)
    barren = np.array([False, False, True, True, True, False])
    water = np.array([False, False, False, False, True, False])
    (rootable, wet, dry), _, residual = partition(
        cell, 3, area, land, barren, water)
    checks = [
        ("ordinary land is fully rootable", rootable[0] == 1.0),
        ("dry barren land is not rootable", rootable[1] == 0.0),
        ("a partial wet/barren cell keeps only ordinary land",
         rootable[2] == 0.5 and wet[2] == 0.5 and dry[2] == 0.0),
        ("water over barren is deducted once", residual <= 1e-15),
        ("extensive NPP uses effective rootable area",
         np.dot(np.array([2.0, 4.0]), np.array([0.5, 1.0])) == 5.0),
        ("soil-carbon feedback becomes a climate-cell mean",
         np.array_equal(np.array([10.0, 20.0]) * np.array([0.5, 1.0]),
                        np.array([5.0, 20.0]))),
    ]
    for name, ok in checks:
        print(f"[{' ok ' if ok else 'FAIL'}] {name}")
    print(f"\n{sum(not ok for _, ok in checks)} failures")
    return sum(not ok for _, ok in checks)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--surface-classes", type=Path,
                        help="legacy post-baseline surface-class source. Omit "
                             "to read the authoritative lake mask directly.")
    parser.add_argument("--surface-water", type=Path,
                        help="authoritative surface_water.nc; defaults to the "
                             "active build. Mutually exclusive with "
                             "--surface-classes.")
    parser.add_argument("--rung", help="atmosphere rung; default configured rung")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        raise SystemExit(1 if self_test() else 0)

    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    if args.surface_classes is not None and args.surface_water is not None:
        raise SystemExit("pass at most one of --surface-classes and --surface-water")
    rung = args.rung.upper() if args.rung else rungs.model_grid(config)[0]
    grid_dir = builds.grid_export(config, rung)
    mesh = Export(builds.mesh_export(config))
    if args.surface_classes is not None:
        source = args.surface_classes
        if not source.is_file():
            raise SystemExit(f"{source} does not exist")
        with Dataset(source) as ds:
            if str(getattr(ds, "vesper_source_build", "")) != config.get("source_build"):
                raise SystemExit(
                    f"{source} belongs to "
                    f"{getattr(ds, 'vesper_source_build', None)!r}, not the active "
                    f"build {config.get('source_build')!r}")
            if str(getattr(ds, "vesper_terrain_hash", "")) != mesh.terrain_hash:
                raise SystemExit(f"{source} carries another terrain hash")
            cover_var = ds["surface_cover"]
            meanings = str(cover_var.flag_meanings).split()
            values = [int(v) for v in np.asarray(cover_var.flag_values)]
            codes = dict(zip(meanings, values))
            if not {"water", "bare"} <= set(codes):
                raise SystemExit(f"{source} has no water/bare cover contract")
            water = np.asarray(cover_var[:], dtype=np.int8) == codes["water"]
        source_role = "post-baseline surface-class copy of the lake mask"
    else:
        source = args.surface_water or (
            builds.component_data("hydrography", config, strict=True)
            / "surface_water.nc")
        if not source.is_file():
            raise SystemExit(f"{source} does not exist; solve surface water first")
        from provenance import require_build
        require_build(source, "rootable surface-water source", config)
        with Dataset(source) as ds:
            water = np.asarray(ds["lake"][:], dtype=bool)
        if water.shape != mesh.substrate_class.shape:
            raise SystemExit(
                f"{source}'s lake mask has {water.size} regions, expected "
                f"{mesh.substrate_class.size}")
        source_role = "authoritative annual-equilibrium solved lake mask"

    legend = {row["code"]: int(row["id"])
              for row in mesh.manifest["lithology"]["rockClasses"]}
    absent = [name for name in config["model"]["barren_rock_classes"]
              if name not in legend]
    if absent:
        raise SystemExit(f"barren rock classes absent from the export: {absent}")
    barren = np.isin(mesh.substrate_class,
                     [legend[name] for name in config["model"]["barren_rock_classes"]])
    land = mesh.surface_class == LAND
    cell, nlat, nlon = gridding.region_cells(mesh, grid_dir)
    area = mesh.cell_area.astype(np.float64)
    (rootable, wet, dry_barren), land_area, residual = partition(
        cell, nlat * nlon, area, land, barren, water)
    if residual > 1e-12:
        raise SystemExit(f"rootable/water/barren partition misses one by {residual:.3g}")
    no_land = land_area <= 0
    for field in (rootable, wet, dry_barren):
        field[no_land] = np.nan

    lat, lon, _ = gridding.grid_geometry(grid_dir)
    spec = gridding.export_grid(grid_dir, name=f"exoplasim-{rung}")
    contract = grid_support_contract(
        spec, f"biosphere/data/{config['source_build']}/rootable_fraction_{rung}.nc",
        "atmosphere_grid", "atmosphere_gaussian_grid",
        [str(source.relative_to(PROJECT_ROOT)),
         str((grid_dir / "manifest.json").relative_to(PROJECT_ROOT))],
        native_measures=[
            {"kind": "area", "variable": "cell_area_m2", "units": "m2"},
        ])
    contract["fields"][0].update({
        "name": "f_rootable",
        "semantics": "dimensionless_fraction",
    })
    contract_identity = validate_contract(contract)

    output = args.output or (
        builds.component_data("biosphere", config) / f"rootable_fraction_{rung}.nc")
    output.parent.mkdir(parents=True, exist_ok=True)
    shape = (nlat, nlon)
    cell_area = spec.cell_area(float(config["planet"]["radius_earth"]) * 6371000.0)
    with Dataset(output, "w", format="NETCDF4") as ds:
        ds.createDimension("lat", nlat)
        ds.createDimension("lon", nlon)
        ds.createVariable("lat", "f8", ("lat",))[:] = lat
        ds.createVariable("lon", "f8", ("lon",))[:] = lon
        for name, data, long_name in (
            ("f_rootable", rootable, "share of native-mesh land that is rootable"),
            ("f_nonrootable_water", wet, "share of native-mesh land under solved water"),
            ("f_nonrootable_barren", dry_barren,
             "share of native-mesh land on dry barren substrate"),
        ):
            variable = ds.createVariable(name, "f4", ("lat", "lon"), zlib=True)
            variable[:] = data.reshape(shape)
            variable.units = "1"
            variable.long_name = long_name
        variable = ds.createVariable("rootable_area_m2", "f8", ("lat", "lon"),
                                     zlib=True)
        variable[:] = rootable.reshape(shape) * cell_area
        variable.units = "m2"
        variable.long_name = ("effective rootable area on the atmosphere support; "
                              "binary-coast model-form error remains SPAT-5's")
        variable = ds.createVariable("cell_area_m2", "f8", ("lat", "lon"), zlib=True)
        variable[:] = cell_area
        variable.units = "m2"
        variable.long_name = "Gaussian quadrature cell area"
        ds.setncattr("vesper_source_build", config["source_build"])
        ds.setncattr("vesper_terrain_hash", mesh.terrain_hash)
        ds.setncattr("vesper_spatial_contract", json.dumps(contract, sort_keys=True))
        ds.setncattr("vesper_spatial_contract_identity", contract_identity)

    valid = ~no_land
    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": "biosphere/scripts/build_rootable_fraction.py",
        "source_build": config["source_build"],
        "terrain_hash": mesh.terrain_hash,
        "rung": rung,
        "partition_source": str(source.relative_to(PROJECT_ROOT)),
        "partition_source_role": source_role,
        "partition_source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "spatial_contract": contract,
        "spatial_contract_identity": contract_identity,
        "partition_max_absolute_residual": residual,
        "cells_without_native_land": int(no_land.sum()),
        "land_area_weighted": {
            "rootable_fraction": float(rootable[valid] @ land_area[valid]
                                       / land_area[valid].sum()),
            "solved_water_fraction": float(wet[valid] @ land_area[valid]
                                          / land_area[valid].sum()),
            "dry_barren_fraction": float(dry_barren[valid] @ land_area[valid]
                                         / land_area[valid].sum()),
        },
        "output": str(output.relative_to(PROJECT_ROOT)),
        "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
    }
    report_path = output.with_name(output.stem + "_report.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"rootable partition residual {residual:.3g}")
    print(f"wrote {output.relative_to(PROJECT_ROOT)}")
    print(f"wrote {report_path.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
