#!/usr/bin/env python3
"""Assess one LPJ-GUESS run before any coupled consumer may read it.

BIO-14: exit code zero is necessary and insufficient. This checks every rank,
output, cell and year; finite and physical values; BIO-12 equilibrium windows;
and cumulative carbon, nitrogen and water closure. The report is written beside
the run manifest and copied into the tracked per-run analysis directory.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import shutil
import struct
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

from _paths import COMPONENT_ROOT, PROJECT_ROOT, RUNS
from lpj_output import (EquilibriumWindowError, reduce_table,
                        timescale_report)

CONFIG = COMPONENT_ROOT / "config" / "lpj_acceptance.yaml"
ANALYSIS = COMPONENT_ROOT / "analysis"
DRIVER_MAGIC = b"VESPDRV8"
DRIVER_HEADER = "<iiiiiidd"
DRIVER_PROVENANCE_BYTES = 64


class AcceptanceError(ValueError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_contract(path: Path = CONFIG) -> dict:
    contract = yaml.safe_load(path.read_text(encoding="utf-8"))
    if contract.get("contract_version") != "vesper-lpj-acceptance/1":
        raise AcceptanceError("unsupported LPJ acceptance contract")
    required = contract.get("required_outputs")
    if not isinstance(required, list) or len(required) != len(set(required)):
        raise AcceptanceError("required_outputs must be a unique list")
    closure = contract.get("closure", {})
    if closure.get("window_complete_forcing_cycles") != 10:
        raise AcceptanceError("closure must use BIO-12's ten complete cycles")
    return contract


def read_table(path: Path) -> tuple[list[str], dict[tuple[float, float, int], np.ndarray]]:
    if not path.is_file():
        raise AcceptanceError(f"missing output {path}")
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise AcceptanceError(f"{path} is empty")
    header = lines[0].split()
    if len(header) < 4 or header[:3] != ["Lon", "Lat", "Year"]:
        raise AcceptanceError(f"{path} must begin Lon Lat Year")
    rows = {}
    for number, line in enumerate(lines[1:], 2):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != len(header):
            raise AcceptanceError(
                f"{path}:{number} has {len(fields)} fields, expected {len(header)}")
        try:
            lon = round(float(fields[0]), 2)
            lat = round(float(fields[1]), 2)
            year_float = float(fields[2])
            year = int(year_float)
            values = np.asarray([float(value) for value in fields[3:]], dtype=float)
        except ValueError as exc:
            raise AcceptanceError(f"{path}:{number} contains nonnumeric data") from exc
        if year_float != year:
            raise AcceptanceError(f"{path}:{number} has nonintegral year")
        if not np.isfinite(values).all():
            raise AcceptanceError(f"{path}:{number} contains NaN or infinity")
        key = (lon, lat, year)
        if key in rows:
            raise AcceptanceError(f"{path}:{number} duplicates {key}")
        rows[key] = values
    if not rows:
        raise AcceptanceError(f"{path} contains no data")
    return header[3:], rows


def coverage(rows: dict, expected_years: int) -> tuple[set[tuple[float, float]], list[int]]:
    cells = {(lon, lat) for lon, lat, _ in rows}
    years = sorted({year for _, _, year in rows})
    if len(years) != expected_years:
        raise AcceptanceError(
            f"output has {len(years)} years, run declares {expected_years}")
    if years != list(range(years[0], years[0] + expected_years)):
        raise AcceptanceError(f"output years are not consecutive: {years}")
    expected = {(lon, lat, year) for lon, lat in cells for year in years}
    missing = expected - set(rows)
    if missing:
        raise AcceptanceError(f"output is missing {len(missing)} cell-year rows")
    return cells, years


def driver_precipitation(path: Path) -> tuple[dict[tuple[float, float], np.ndarray], dict]:
    sidecar = path.with_name(path.stem + "_provenance.json")
    if not sidecar.is_file():
        raise AcceptanceError(f"{path} has no provenance sidecar")
    provenance = json.loads(sidecar.read_text(encoding="utf-8"))
    layers = provenance.get("physical_layers")
    if not isinstance(layers, int) or layers < 1:
        raise AcceptanceError(f"{sidecar} has no physical_layers")
    data = path.read_bytes()
    header_size = struct.calcsize(DRIVER_HEADER)
    if data[:8] != DRIVER_MAGIC or len(data) < 8 + header_size:
        raise AcceptanceError(f"{path} is not a complete VESPDRV8 driver")
    ncells, nintervals, year_days, nyears, subdaily, padding, co2, ndep = \
        struct.unpack(DRIVER_HEADER, data[8:8 + header_size])
    if subdaily != 0 or padding != 0:
        raise AcceptanceError("BIO-14 cannot parse a sampled or unaligned driver")
    offset = 8 + header_size + DRIVER_PROVENANCE_BYTES
    interval_count = nyears * nintervals * 4
    intervals = np.frombuffer(data, dtype="<f8", count=interval_count,
                              offset=offset).reshape(nyears, nintervals, 4)
    durations_days = intervals[:, :, 2] / 86400.0
    if not np.all(durations_days > 0.0):
        raise AcceptanceError("driver has a non-positive interval duration")
    offset += interval_count * 8
    nvalues = nyears * nintervals
    record_bytes = 16 + 8 + (4 + layers) * 8 + 4 * nvalues * 8
    result = {}
    for _ in range(ncells):
        if offset + record_bytes > len(data):
            raise AcceptanceError(f"{path} ends inside a cell record")
        lon, lat = struct.unpack_from("<dd", data, offset)
        offset += 16 + 8 + (4 + layers) * 8
        offset += nvalues * 8  # temperature
        precip = np.frombuffer(data, dtype="<f8", count=nvalues,
                               offset=offset).reshape(nyears, nintervals).copy()
        offset += nvalues * 8
        offset += 2 * nvalues * 8  # shortwave and temperature range
        result[(round(lon, 2), round(lat, 2))] = np.sum(
            precip * durations_days, axis=1)
    if offset != len(data):
        raise AcceptanceError(f"{path} has {len(data) - offset} trailing bytes")
    return result, {"cells": ncells, "cycle_years": nyears,
                    "year_length_days": year_days, "co2_ppm": co2,
                    "ndep_kg_n_ha_earth_year": ndep,
                    "sha256": sha256(path), "provenance_sha256": sha256(sidecar)}


def check_physical(name: str, fields: list[str], rows: dict, contract: dict) -> None:
    values = np.stack(list(rows.values()))
    if name in contract.get("nonnegative_outputs", []):
        tolerance = float(contract["negative_roundoff_tolerance"])
        if float(values.min()) < -tolerance:
            raise AcceptanceError(f"{name} has a negative physical value {values.min()}")
    bound = contract.get("fraction_outputs", {}).get(name)
    if isinstance(bound, list):
        low, high = map(float, bound)
        if float(values.min()) < low or float(values.max()) > high:
            raise AcceptanceError(f"{name} leaves [{low}, {high}]")
    elif isinstance(bound, dict):
        for field, limits in bound.items():
            if field not in fields:
                raise AcceptanceError(f"{name} has no bounded field {field}")
            column = values[:, fields.index(field)]
            low, high = map(float, limits)
            if float(column.min()) < low or float(column.max()) > high:
                raise AcceptanceError(f"{name}.{field} leaves [{low}, {high}]")


def field_series(tables: dict, output: str, field: str,
                 cell: tuple[float, float], years: list[int]) -> np.ndarray:
    fields, rows = tables[output]
    if field not in fields:
        raise AcceptanceError(f"{output} has no {field} field")
    index = fields.index(field)
    return np.asarray([rows[(cell[0], cell[1], year)][index] for year in years])


def closure_report(tables: dict, cells: set, years: list[int], manifest: dict,
                   contract: dict) -> dict:
    closure = contract["closure"]
    cycles = int(closure["window_complete_forcing_cycles"])
    cycle_years = int(manifest.get("forcing", {}).get("cycle_years", 1))
    count = cycles * cycle_years
    if len(years) < count:
        raise AcceptanceError(f"closure needs {count} end years")
    selected = years[-count:]
    driver_record = manifest.get("inputs", {}).get("driver", {})
    driver = Path(driver_record.get("path", ""))
    if not driver.is_file() or sha256(driver) != driver_record.get("sha256"):
        raise AcceptanceError("manifest-pinned driver is absent or changed")
    precip, precip_identity = driver_precipitation(driver)
    if set(precip) != cells:
        raise AcceptanceError(
            f"driver/output support differs: {len(precip)} versus {len(cells)} cells")

    reports = {}
    for element in ("carbon", "nitrogen"):
        rule = closure[element]
        residuals = []
        throughputs = []
        for cell in cells:
            pool = field_series(tables, rule["pool_output"], rule["pool_field"],
                                cell, selected)
            flux = field_series(tables, rule["flux_output"], rule["flux_field"],
                                cell, selected)
            factor = float(rule.get("pool_to_flux_units", 1.0))
            residual = (pool[-1] - pool[0]) * factor + float(flux[1:].sum())
            throughput = float(np.abs(flux[1:]).sum())
            limit = max(float(rule[next(k for k in rule if k.startswith("absolute_floor_") )]),
                        float(rule["relative_throughput_limit"]) * throughput)
            residuals.append(residual)
            throughputs.append((throughput, limit))
        failed = [abs(value) > limit for value, (_, limit)
                  in zip(residuals, throughputs)]
        reports[element] = {
            "maximum_absolute_residual": float(np.max(np.abs(residuals))),
            "maximum_allowed_residual": float(max(limit for _, limit in throughputs)),
            "failed_cells": int(sum(failed)), "cells": len(cells),
        }
        if any(failed):
            raise AcceptanceError(f"{element} closure fails in {sum(failed)} cells")

    water = closure["water"]
    residuals = []
    limits = []
    for cell in cells:
        aet = field_series(tables, water["aet_output"], water["aet_field"],
                           cell, selected)
        runoff = field_series(tables, water["runoff_output"], water["runoff_field"],
                              cell, selected)
        p_cycle = float(precip[cell].sum())
        p_total = cycles * p_cycle
        residual = p_total - float(aet.sum()) - float(runoff.sum())
        limit = max(float(water["absolute_floor_mm"]),
                    float(water["relative_throughput_limit"]) * p_total)
        residuals.append(residual)
        limits.append(limit)
    failed = [abs(value) > limit for value, limit in zip(residuals, limits)]
    reports["water"] = {
        "form": "precipitation - AET - runoff = implied end-minus-start storage",
        "interpretation": water["interpretation"],
        "maximum_absolute_residual_mm": float(np.max(np.abs(residuals))),
        "maximum_allowed_residual_mm": float(max(limits)),
        "failed_cells": int(sum(failed)), "cells": len(cells),
        "driver": precip_identity,
    }
    if any(failed):
        raise AcceptanceError(f"water closure fails in {sum(failed)} cells")
    return reports


def assess(run_dir: Path, *, contract_path: Path = CONFIG,
           write: bool = True) -> dict:
    run_dir = Path(run_dir).resolve()
    contract = read_contract(contract_path)
    manifest_path = run_dir / "run_manifest.json"
    if not manifest_path.is_file():
        raise AcceptanceError(f"{run_dir} has no run_manifest.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    physical = manifest.get("physical", {})
    ranks = int(physical.get("ranks", manifest.get("ranks", 0)))
    nyear = int(physical.get("nyear", 0))
    if ranks < 1 or nyear < 1:
        raise AcceptanceError("manifest has invalid ranks or nyear")

    tables = {}
    canonical_support = None
    canonical_years = None
    output_hashes = {}
    rank_summary = {}
    for output in contract["required_outputs"]:
        merged_fields, merged_rows = read_table(run_dir / output)
        cells, years = coverage(merged_rows, nyear)
        check_physical(output, merged_fields, merged_rows, contract)
        pieces = []
        rank_cells = []
        for rank in range(1, ranks + 1):
            fields, rows = read_table(run_dir / f"run{rank}" / output)
            if fields != merged_fields:
                raise AcceptanceError(f"rank {rank} {output} has a different header")
            piece_cells, piece_years = coverage(rows, nyear)
            if piece_years != years:
                raise AcceptanceError(f"rank {rank} {output} has different years")
            pieces.append(rows)
            rank_cells.append(piece_cells)
        for left in range(ranks):
            for right in range(left + 1, ranks):
                if rank_cells[left] & rank_cells[right]:
                    raise AcceptanceError(f"{output} repeats cells across ranks")
        union = {}
        for piece in pieces:
            union.update(piece)
        if set(union) != set(merged_rows):
            raise AcceptanceError(f"merged {output} is not the exact rank union")
        for key in union:
            if not np.array_equal(union[key], merged_rows[key]):
                raise AcceptanceError(f"merged {output} changes rank row {key}")
        counts = [len(value) for value in rank_cells]
        if max(counts) - min(counts) > 1:
            raise AcceptanceError(f"{output} rank cell counts are imbalanced: {counts}")
        if canonical_support is None:
            canonical_support, canonical_years = cells, years
        elif cells != canonical_support or years != canonical_years:
            raise AcceptanceError(f"{output} has different cell/year support")
        tables[output] = (merged_fields, merged_rows)
        output_hashes[output] = sha256(run_dir / output)
        rank_summary[output] = counts

    forcing_cells = int(manifest.get("forcing", {}).get("cells", 0))
    if forcing_cells != len(canonical_support):
        raise AcceptanceError(
            f"forcing declares {forcing_cells} cells, outputs have {len(canonical_support)}")

    stability = {}
    for output in contract["stability_outputs"]:
        try:
            reduced = reduce_table(run_dir / output)
        except EquilibriumWindowError as exc:
            raise AcceptanceError(f"{output} equilibrium refusal: {exc}") from exc
        stability[output] = {
            "window": reduced.report.get("window"),
            "trend": reduced.report.get("trend"),
        }
    timescales = _timescales_or_reason(run_dir)
    closures = closure_report(tables, canonical_support, canonical_years,
                              manifest, contract)
    report = {
        "contract_version": contract["contract_version"],
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "run_id": manifest.get("run_id", run_dir.name),
        "source_build": manifest.get("source_build"),
        "manifest_sha256": sha256(manifest_path),
        "contract_sha256": sha256(contract_path),
        "timescales": timescales,
        "coverage": {"ranks": ranks, "cells": len(canonical_support),
                     "years": canonical_years,
                     "rank_cell_counts_by_output": rank_summary},
        "output_sha256": output_hashes,
        "stability": stability,
        "closure": closures,
        "verdict": "PASS",
    }
    if write:
        target = run_dir / "acceptance.json"
        target.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        analysis_dir = ANALYSIS / report["run_id"]
        analysis_dir.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(target, analysis_dir / "acceptance.json")
    return report


def _timescales_or_reason(run_dir: Path) -> dict:
    """The ecological timescales, or why they could not be taken."""
    try:
        return timescale_report(run_dir)
    except (EquilibriumWindowError, OSError, ValueError, KeyError) as exc:
        return {"measured": False, "reason": str(exc)}


def write_failure(run_dir: Path, error: Exception, *,
                  contract_path: Path = CONFIG) -> dict:
    """Persist a refusal even when assessment cannot build a PASS report."""
    run_dir = Path(run_dir).resolve()
    manifest_path = run_dir / "run_manifest.json"
    manifest = {}
    if manifest_path.is_file():
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            pass
    report = {
        "contract_version": "vesper-lpj-acceptance/1",
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "run_id": manifest.get("run_id", run_dir.name),
        "source_build": manifest.get("source_build"),
        "manifest_sha256": sha256(manifest_path) if manifest_path.is_file() else None,
        "contract_sha256": sha256(contract_path) if contract_path.is_file() else None,
        "verdict": "FAIL",
        "refusal": str(error),
        # A run refused for not having settled is exactly the run whose
        # timescales say how long the next one has to be, so they are recorded on
        # the refusal and not only on a pass. `lib/run_lengths.py` reads them from
        # here. Guarded, because a refusal raised before the tables were readable
        # must not be masked by a second failure while measuring them.
        "timescales": _timescales_or_reason(run_dir),
    }
    run_dir.mkdir(parents=True, exist_ok=True)
    target = run_dir / "acceptance.json"
    target.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    analysis_dir = ANALYSIS / report["run_id"]
    analysis_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(target, analysis_dir / "acceptance.json")
    return report


def _write_driver(path: Path, cells: list[tuple[float, float]], precip_mm: float) -> None:
    nyears, nintervals, year_days, layers = 1, 1, 183, 2
    interval = np.asarray([[[0.0, year_days * 86400.0,
                             year_days * 86400.0, 0.0]]], dtype="<f8")
    with path.open("wb") as handle:
        handle.write(DRIVER_MAGIC)
        handle.write(struct.pack(DRIVER_HEADER, len(cells), nintervals, year_days,
                                 nyears, 0, 0, 450.0, 1.0))
        handle.write(b"0" * DRIVER_PROVENANCE_BYTES)
        handle.write(interval.tobytes())
        for lon, lat in cells:
            handle.write(struct.pack("<ddii", lon, lat, 2, 0))
            handle.write(struct.pack("<6d", 1.0, 0.4, 0.3, 0.1, 1.0, 1.0))
            rate = precip_mm / year_days
            for value in (280.0, rate, 100.0, 5.0):
                handle.write(struct.pack("<d", value))
    path.with_name(path.stem + "_provenance.json").write_text(
        json.dumps({"physical_layers": layers}) + "\n", encoding="utf-8")


def _table_text(output: str, cells: list[tuple[float, float]], years: range) -> str:
    special = {
        "cpool.out": (["Total"], [10.0]), "cflux.out": (["NEE"], [0.0]),
        "npool.out": (["Total"], [0.1]), "nflux.out": (["NEE"], [0.0]),
        "aaet.out": (["Total"], [600.0]),
        "tot_runoff.out": (["Surf", "Drain", "Base", "Total"],
                            [100.0, 200.0, 100.0, 400.0]),
        "fpc.out": (["Tree", "Total"], [0.4, 0.4]),
        "firert.out": (["FireRT", "BurntFr"], [100.0, 0.01]),
        "ngases.out": (["NH3_soil", "Total"], [0.01, 0.01]),
    }
    fields, values = special.get(output, (["Total"], [1.0]))
    lines = ["Lon Lat Year " + " ".join(fields)]
    for lon, lat in cells:
        for year in years:
            lines.append(f"{lon} {lat} {year} " + " ".join(map(str, values)))
    return "\n".join(lines) + "\n"


def selftest() -> dict:
    fixtures = []
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        run = root / "run"
        run.mkdir()
        cells = [(0.0, 10.0), (20.0, -10.0)]
        years = range(12)
        driver = root / "driver.bin"
        _write_driver(driver, cells, 1000.0)
        contract = read_contract()
        for rank in (1, 2):
            (run / f"run{rank}").mkdir()
        for output in contract["required_outputs"]:
            parts = []
            for rank, cell in enumerate(cells, 1):
                text = _table_text(output, [cell], years)
                (run / f"run{rank}" / output).write_text(text, encoding="utf-8")
                parts.append(text.splitlines())
            merged = parts[0] + parts[1][1:]
            (run / output).write_text("\n".join(merged) + "\n", encoding="utf-8")
        manifest = {
            "run_id": "fixture", "source_build": "fixture",
            "physical": {"nyear": 12, "ranks": 2}, "ranks": 2,
            "forcing": {"cells": 2, "cycle_years": 1},
            "inputs": {"driver": {"path": str(driver), "sha256": sha256(driver)}},
        }
        (run / "run_manifest.json").write_text(json.dumps(manifest) + "\n")
        assess(run, write=False)
        fixtures.append({"fixture": "complete run", "pass": True})

        mutations = {
            "missing rank output": lambda bed: (bed / "run2" / "lai.out").unlink(),
            "non-finite value": lambda bed: (bed / "anpp.out").write_text(
                (bed / "anpp.out").read_text().replace("1.0", "nan", 1)),
            "negative stock": lambda bed: (bed / "cpool.out").write_text(
                (bed / "cpool.out").read_text().replace("10.0", "-1.0", 1)),
            "carbon leak": lambda bed: (bed / "cflux.out").write_text(
                (bed / "cflux.out").read_text().replace("0.0", "1.0")),
            "water leak": lambda bed: (bed / "aaet.out").write_text(
                (bed / "aaet.out").read_text().replace("600.0", "500.0")),
        }
        for name, mutate in mutations.items():
            bed = root / name.replace(" ", "_")
            shutil.copytree(run, bed)
            mutate(bed)
            try:
                assess(bed, write=False)
            except AcceptanceError:
                fixtures.append({"fixture": name, "pass": True})
            else:
                fixtures.append({"fixture": name, "pass": False})
    return {"fixtures": fixtures,
            "verdict": "PASS" if all(row["pass"] for row in fixtures) else "FAIL"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", nargs="?", type=Path)
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()
    if args.selftest:
        report = selftest()
        for row in report["fixtures"]:
            print(f"[{' ok ' if row['pass'] else 'FAIL'}] {row['fixture']}")
        raise SystemExit(0 if report["verdict"] == "PASS" else 1)
    if args.run is None:
        raise SystemExit("provide a run directory or --selftest")
    try:
        report = assess(args.run)
    except (AcceptanceError, EquilibriumWindowError, OSError,
            json.JSONDecodeError, yaml.YAMLError) as exc:
        report = write_failure(args.run, exc)
        raise SystemExit(f"FAIL: {report['run_id']}: {report['refusal']}") from exc
    print(f"PASS: {report['run_id']} ({report['coverage']['cells']} cells, "
          f"{len(report['coverage']['years'])} years, "
          f"{report['coverage']['ranks']} ranks)")


if __name__ == "__main__":
    main()
