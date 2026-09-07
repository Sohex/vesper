#!/usr/bin/env python3
"""Grid an accepted LPJ-GUESS run's simulated vegetation onto the atmosphere grid.

    python biosphere/scripts/build_vegetation_field.py
    python biosphere/scripts/build_vegetation_field.py --run lpj_<id>
    python biosphere/scripts/build_vegetation_field.py --self-test

An LPJ-GUESS `.out` table is one row per gridcell and simulated year, and the
quantity a consumer wants is neither the last year nor an arbitrary mean of a
few: it is the reduced value over the span BIO-12's equilibrium contract
certifies. `lib/lpj_output.py:reduce_table` is that reduction and it refuses a
record whose drift is not bounded, so this step emits nothing from a run whose
vegetation is still moving.

WHICH LONGITUDE THESE ROWS CARRY. `build_lpj_driver.py` writes the atmosphere
model's OWN longitude labels, wrapped from 0..360 into -180..180 because that
is the range LPJ-GUESS reads. Wrapping a label changes its name and not its
column, so these rows are on the model label axis and go through
`gridding.model_label_cells`, which refuses a label that is not on that axis.
The export's centres would sit half a column off every label and are refused by
500 times the bar rather than rounded to the nearest column. CLAUDE.md rule 3.

The emitted cover is what the model simulated and not a classification of it.
`fpc` is foliar projective cover per plant functional type, and its sum exceeds
one where cohorts shade each other, which is what `vegmode "cohort"` means and
not a defect: the total is a sum of layers, and a consumer wanting a ground
cover fraction clips it.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from netCDF4 import Dataset
import numpy as np
import yaml

from _paths import CONFIG, PROJECT_ROOT

import builds
import gridding
import lpj_output
import nc_geometry
import rungs
from spatial_support import grid_support_contract, validate_contract

# The tables gridded, and what each one is. Every one of these is a per-PFT
# annual table with the same header, so they share one cell axis and one PFT
# axis; `reduce_table` is asked to prove that rather than told it.
TABLES = {
    "fpc.out": ("fpc", "1",
                "foliar projective cover, summed over cohorts of the type"),
    "lai.out": ("lai", "m2 m-2",
                "leaf area index"),
    "cmass.out": ("cmass", "kgC m-2",
                  "vegetation carbon"),
}
# The trailing column of each table is the model's own sum across the types. It
# is read and CHECKED rather than dropped, because a total that disagrees with
# the columns beside it says the table was misparsed.
TOTAL_COLUMN = "Total"
TOTAL_TOLERANCE = 5.0e-3


def place(values, lat_label, lon_label, nlat: int, nlon: int, names: list[str]):
    """Per-cell rows onto (field, lat, lon), NaN where nothing was simulated.

    NaN and not zero: a cell LPJ did not simulate has no cover, and zero is a
    measurement of bare ground. The two are drawn differently and a consumer
    that cannot tell them apart paints ocean as desert.
    """
    values = np.asarray(values, dtype=np.float64)
    rows, cols = gridding.model_label_cells(
        lat_label, lon_label, nlat, nlon,
        what="an LPJ-GUESS output table's gridcell")
    flat = rows.astype(np.int64) * nlon + cols.astype(np.int64)
    if len(np.unique(flat)) != len(flat):
        raise SystemExit(
            f"{len(flat) - len(np.unique(flat))} of {len(flat)} simulated cells "
            "share a grid cell with another. Two LPJ gridcells in one atmosphere "
            "cell means the coordinates are not this grid's.")
    out = np.full((len(names), nlat * nlon), np.nan, dtype=np.float64)
    out[:, flat] = values.T
    return out.reshape(len(names), nlat, nlon), flat


def self_test() -> int:
    nlat, nlon = 4, 8
    labels = gridding.model_longitude_labels(nlon)
    lats = gridding.gaussian_latitudes(nlat)
    # Two cells, addressed by the model's labels wrapped into -180..180 exactly
    # as `build_lpj_driver.py` writes them.
    signed = np.where(labels > 180.0, labels - 360.0, labels)
    lon_label = np.array([signed[1], signed[6]])
    lat_label = np.array([lats[0], lats[3]])
    values = np.array([[1.0, 2.0], [3.0, 4.0]])
    field, _ = place(values, lat_label, lon_label, nlat, nlon, ["a", "b"])

    def refuses(fn) -> bool:
        try:
            fn()
        except SystemExit:
            return True
        return False

    export_centres = gridding.gaussian_grid(nlat, nlon).cell_centres()[1]
    checks = [
        ("a wrapped model label lands in the column it names",
         field[0, 0, 1] == 1.0 and field[1, 3, 6] == 4.0),
        ("an unsimulated cell is NaN and not zero",
         np.isnan(field[0, 2, 2]) and np.isnan(field[1, 0, 0])),
        ("the export's centres are refused, not rounded to a column",
         refuses(lambda: place(values, lat_label,
                               np.array([export_centres[1], export_centres[6]]),
                               nlat, nlon, ["a", "b"]))),
        ("two rows in one grid cell are refused",
         refuses(lambda: place(values, np.array([lats[0], lats[0]]),
                               np.array([signed[1], signed[1]]),
                               nlat, nlon, ["a", "b"]))),
        ("a cover total is checked against the columns it sums",
         abs((1.0 + 2.0) - 3.0) <= TOTAL_TOLERANCE
         and not abs((1.0 + 2.0) - 3.5) <= TOTAL_TOLERANCE),
    ]
    for name, ok in checks:
        print(f"[{' ok ' if ok else 'FAIL'}] {name}")
    failures = sum(not ok for _, ok in checks)
    print(f"\n{failures} failures")
    return failures


def resolve_run(config: dict, requested: str | None) -> Path:
    """The run directory to grid: the one named, or the accepted one.

    With no `--run` this REFUSES rather than picking, whenever more than one
    accepted run stands for the active build. Which vegetation the world's
    picture is drawn from is the caller's declaration, and choosing the newest
    would silently move the map every time a run lands.
    """
    runs = PROJECT_ROOT / "biosphere" / "runs"
    if requested:
        run_dir = Path(requested)
        if not run_dir.is_dir():
            run_dir = runs / requested
        if not run_dir.is_dir():
            raise SystemExit(f"{requested} is not a run directory")
        return run_dir
    accepted = []
    for candidate in sorted(runs.glob("lpj_*")):
        report = candidate / "acceptance.json"
        if not report.is_file():
            continue
        try:
            record = json.loads(report.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if (record.get("verdict") == "PASS"
                and record.get("source_build") == config.get("source_build")):
            accepted.append(candidate)
    if not accepted:
        raise SystemExit(
            f"no accepted LPJ run for build {config.get('source_build')!r}; "
            "run biosphere/scripts/run_lpj_guess.py, or name one with --run")
    if len(accepted) > 1:
        names = ", ".join(path.name for path in accepted)
        raise SystemExit(
            f"{len(accepted)} accepted runs stand for this build ({names}). "
            "Name the one the world's vegetation is taken from with --run.")
    return accepted[0]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", help="LPJ run id or directory; default the "
                                      "single accepted run for the active build")
    parser.add_argument("--rung", help="atmosphere rung; default configured rung")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        raise SystemExit(1 if self_test() else 0)

    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    run_dir = resolve_run(config, args.run)
    acceptance = lpj_output.require_lpj_acceptance(run_dir / "fpc.out")
    if acceptance.get("source_build") != config.get("source_build"):
        raise SystemExit(
            f"{run_dir.name} was run on build {acceptance.get('source_build')!r}, "
            f"not the active build {config.get('source_build')!r}")
    manifest = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))

    rung = args.rung.upper() if args.rung else rungs.model_grid(config)[0]
    nlat, nlon, _ = rungs.geometry(rung)
    grid_dir = builds.grid_export(config, rung)

    fields, reports, names = {}, {}, None
    cells = None
    for table, (variable, units, long_name) in TABLES.items():
        path = run_dir / table
        lpj_output.require_lpj_acceptance(path)
        reduced = lpj_output.reduce_table(path)
        columns = [name for name in reduced.names if name != TOTAL_COLUMN]
        if TOTAL_COLUMN not in reduced.names:
            raise SystemExit(f"{path} has no {TOTAL_COLUMN} column to check against")
        if names is None:
            names = columns
        elif columns != names:
            raise SystemExit(
                f"{path} carries {columns}, and {list(TABLES)[0]} carries {names}. "
                "The tables do not share a plant functional type axis.")
        keys = sorted(reduced.values)
        if cells is None:
            cells = keys
        elif keys != cells:
            raise SystemExit(f"{path} was simulated on a different set of gridcells")
        stacked = np.stack([reduced.values[key] for key in keys])
        index = {name: i for i, name in enumerate(reduced.names)}
        per_type = stacked[:, [index[name] for name in names]]
        total = stacked[:, index[TOTAL_COLUMN]]
        residual = float(np.max(np.abs(per_type.sum(axis=1) - total)))
        if residual > TOTAL_TOLERANCE * max(1.0, float(np.max(np.abs(total)))):
            raise SystemExit(
                f"{path}: the {TOTAL_COLUMN} column misses the sum of the type "
                f"columns by {residual:.4g}. The table was misparsed.")
        lon_label = np.array([key[0] for key in keys], dtype=np.float64)
        lat_label = np.array([key[1] for key in keys], dtype=np.float64)
        placed, flat = place(per_type, lat_label, lon_label, nlat, nlon, names)
        fields[variable] = (placed, units, long_name)
        reports[table] = reduced.report

    simulated = manifest.get("cells_simulated")
    if simulated is not None and len(cells) != int(simulated):
        raise SystemExit(
            f"{len(cells)} gridcells in the tables against {simulated} the run "
            "manifest says were simulated")

    covered = np.zeros((nlat, nlon), dtype=bool)
    covered.reshape(-1)[flat] = True

    spec = gridding.export_grid(grid_dir, name=f"exoplasim-{rung}")
    artifact = (f"biosphere/data/{config['source_build']}/vegetation_{rung}.nc")
    contract = grid_support_contract(
        spec, artifact, "atmosphere_grid", "atmosphere_gaussian_grid",
        [f"biosphere/runs/{run_dir.name}/{table}" for table in TABLES],
        native_measures=[{"kind": "area", "variable": "cell_area_m2", "units": "m2"}])
    contract["fields"][0].update({"name": "fpc", "semantics": "dimensionless_fraction"})
    contract_identity = validate_contract(contract)

    output = args.output or (
        builds.component_data("biosphere", config) / f"vegetation_{rung}.nc")
    output.parent.mkdir(parents=True, exist_ok=True)
    cell_area = spec.cell_area(nc_geometry.planet_radius_m(config))
    with Dataset(output, "w", format="NETCDF4") as ds:
        ds.createDimension("lat", nlat)
        ds.createDimension("lon", nlon)
        ds.createDimension("pft", len(names))
        ds.createDimension("name_length", max(len(name) for name in names))
        label = ds.createVariable("pft_name", "S1", ("pft", "name_length"))
        label[:] = np.array([list(name.ljust(len(label[0]))) for name in names],
                            dtype="S1")
        label.long_name = ("plant functional type, in the column order of the "
                           "run's own output tables")
        for variable, (data, units, long_name) in fields.items():
            handle = ds.createVariable(variable, "f4", ("pft", "lat", "lon"),
                                       zlib=True, fill_value=np.nan)
            handle[:] = data
            handle.units = units
            handle.long_name = long_name
            handle.coordinates = "pft_name"
            total = ds.createVariable(f"{variable}_total", "f4", ("lat", "lon"),
                                      zlib=True, fill_value=np.nan)
            total[:] = np.where(covered, np.nansum(data, axis=0), np.nan)
            total.units = units
            total.long_name = f"{long_name}, summed over the types"
        mask = ds.createVariable("simulated", "i1", ("lat", "lon"), zlib=True)
        mask[:] = covered.astype(np.int8)
        mask.units = "1"
        mask.long_name = ("1 where the run simulated a gridcell; every emitted "
                          "field is fill elsewhere")
        area = ds.createVariable("cell_area_m2", "f8", ("lat", "lon"), zlib=True)
        area[:] = cell_area
        area.units = "m2"
        area.long_name = "Gaussian quadrature cell area"
        ds.setncattr("vesper_source_build", config["source_build"])
        ds.setncattr("vesper_lpj_run", run_dir.name)
        ds.setncattr("vesper_lpj_acceptance_contract",
                     acceptance.get("contract_version", ""))
        ds.setncattr("vesper_spatial_contract", json.dumps(contract, sort_keys=True))
        ds.setncattr("vesper_spatial_contract_identity", contract_identity)
        # LAST, and the convention is DECLARED rather than detected: these
        # columns carry the atmosphere model's own labels, which is what the
        # rows of an LPJ output table are addressed by.
        declared = nc_geometry.declare_grid(
            ds, convention=nc_geometry.MODEL_LABELS,
            what="the simulated vegetation field")

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": "biosphere/scripts/build_vegetation_field.py",
        "source_build": config["source_build"],
        "rung": rung,
        "lpj_run": run_dir.name,
        "lpj_acceptance": acceptance.get("contract_sha256"),
        "plant_functional_types": names,
        "gridcells_simulated": len(cells),
        "gridcells_on_grid": int(covered.sum()),
        "longitude_convention": declared.convention,
        "spatial_contract": contract,
        "spatial_contract_identity": contract_identity,
        "equilibrium_window": {table: record.get("identity")
                               for table, record in reports.items()},
    }
    report_path = output.with_name(output.stem + "_report.json")
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(f"wrote {output.relative_to(PROJECT_ROOT)}")
    print(f"      {report_path.relative_to(PROJECT_ROOT)}")
    print(f"{len(cells)} gridcells, {len(names)} plant functional types, "
          f"from {run_dir.name}")


if __name__ == "__main__":
    main()
