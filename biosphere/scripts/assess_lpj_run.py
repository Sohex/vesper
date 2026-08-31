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

# `_paths` is what puts lib/ on the path, so it is imported before anything
# that lives there.
from _paths import COMPONENT_ROOT, PROJECT_ROOT, RUNS

import lpj_table
from lpj_output import (EquilibriumWindowError, reduce_table,
                        timescale_report)

CONFIG = COMPONENT_ROOT / "config" / "lpj_acceptance.yaml"
# The equilibrium contract, read HERE only so the self-test's fixture tables
# carry the columns it names. `lib/lpj_output.py` is what applies it.
EQUILIBRIUM_CONFIG = COMPONENT_ROOT / "config" / "equilibrium_window.yaml"
ANALYSIS = COMPONENT_ROOT / "analysis"
DRIVER_MAGIC = b"VESPDRV8"
DRIVER_HEADER = "<iiiiiidd"
DRIVER_PROVENANCE_BYTES = 64


class AcceptanceError(ValueError):
    pass


class ClosureError(AcceptanceError):
    """A closure refusal that carries the per-cell report it refused on.

    A refusal that says four cells failed without saying WHICH is a refusal
    nobody can act on, and the closure check raises before `assess` has built
    anything to write. So the report travels on the exception and
    `write_failure` records it beside the refusal text.
    """

    def __init__(self, message: str, closure: dict):
        super().__init__(message)
        self.closure = closure


# How many failing cells a refusal names. A refusal has to name enough of them
# to be actionable and not so many that the report becomes the output file
# again; the ones it names are the largest residuals, and the count is always
# reported in full.
NAMED_FAILURES = 64
# The most decimal places a written column is searched for. Past this a column
# is treated as unquantised, which is the benign direction: it makes the
# resolution bound smaller, never larger.
MAX_WRITTEN_DECIMALS = 12


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


def read_table(path: Path) -> lpj_table.Table:
    """One output table, read column-wise wherever that can be certified.

    `lib/lpj_table.py` parses each column in one pass and declines anything it
    cannot certify; the row-at-a-time parser below is what then diagnoses the
    defect, so every refusal a malformed output earns is worded exactly as
    BIO-14 has always worded it.
    """
    if not path.is_file():
        raise AcceptanceError(f"missing output {path}")
    try:
        return lpj_table.read(path)
    except lpj_table.RowParseRequired:
        return read_table_rows(path)


def read_table_rows(path: Path) -> lpj_table.Table:
    """The row-at-a-time parse: the only thing that diagnoses a bad table."""
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
    return _table_from_rows(path, header[3:], rows)


def _table_from_rows(path: Path, fields: list[str], rows: dict) -> lpj_table.Table:
    """A row-parsed table in the sorted, integer-keyed form every check reads."""
    ordered = sorted(rows)
    lon = np.asarray([key[0] for key in ordered], dtype=float)
    lat = np.asarray([key[1] for key in ordered], dtype=float)
    year = np.asarray([key[2] for key in ordered], dtype=np.int64)
    if year.size and np.any(np.abs(year) >= lpj_table.YEAR_OFFSET):
        raise AcceptanceError(f"{path} has a year outside the assessable range")
    key = lpj_table.composite_key(np.rint(lon * lpj_table.CENTS).astype(np.int64),
                                  np.rint(lat * lpj_table.CENTS).astype(np.int64),
                                  year)
    return lpj_table.Table(fields=fields, lon=lon, lat=lat, year=year,
                           values=np.stack([rows[name] for name in ordered]),
                           key=key)


def coverage(table: lpj_table.Table,
             expected_years: int) -> tuple[list[tuple[float, float]], list[int]]:
    """The cells and years an output covers, refusing an incomplete product.

    A table's rows are unique by construction, so the product is complete
    exactly when the row count is the product of the two supports, and the
    shortfall is how many cell-year rows are missing.
    """
    years = table.years()
    if len(years) != expected_years:
        raise AcceptanceError(
            f"output has {len(years)} years, run declares {expected_years}")
    if years != list(range(years[0], years[0] + expected_years)):
        raise AcceptanceError(f"output years are not consecutive: {years}")
    cells = table.cells()
    missing = len(cells) * len(years) - table.rows
    if missing:
        raise AcceptanceError(f"output is missing {missing} cell-year rows")
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


def check_physical(name: str, fields: list[str], values: np.ndarray,
                   contract: dict) -> None:
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


def verify_merge(output: str, merged: lpj_table.Table,
                 pieces: list[lpj_table.Table]) -> None:
    """The merged table is the rank union, cell-year for cell-year and value for value.

    Every table arrives sorted on one integer cell-year key, so the union is a
    concatenation put back into that order and the whole check is two array
    comparisons. The ranks have already been shown to hold disjoint cells, so
    no key appears twice in the concatenation and comparing the sorted key
    arrays is the same statement as comparing the two supports as sets.
    """
    keys = np.concatenate([piece.key for piece in pieces])
    order = np.argsort(keys, kind="stable")
    if not np.array_equal(keys[order], merged.key):
        raise AcceptanceError(f"merged {output} is not the exact rank union")
    values = np.concatenate([piece.values for piece in pieces])[order]
    if np.array_equal(values, merged.values):
        return
    row = int(np.flatnonzero(np.any(values != merged.values, axis=1))[0])
    key = (float(merged.lon[row]), float(merged.lat[row]), int(merged.year[row]))
    raise AcceptanceError(f"merged {output} changes rank row {key}")


def field_series(tables: dict, output: str, field: str, cell: int,
                 years: range) -> np.ndarray:
    """One cell's series for one field, over a contiguous run of year indices.

    The cell and the years are INDEX positions into the cell-year grid the
    coverage check has already established is complete, so the series is read
    straight out of that grid rather than looked up a year at a time.
    """
    fields, grid = tables[output]
    if field not in fields:
        raise AcceptanceError(f"{output} has no {field} field")
    index = fields.index(field)
    return np.asarray([grid[cell, year, index] for year in years])


def written_quantum(tables: dict, output: str, field: str) -> float:
    """The decimal step one column of an output table is WRITTEN on.

    LPJ-GUESS writes each column at a compiled-in fixed precision
    (`modules/commonoutput.cpp`, `ColumnDescriptor(width, precision)`), so a
    value read back from the table is the simulated quantity rounded to that
    many decimals and the column cannot resolve anything finer. This measures
    that step from the ARTIFACT rather than from the model source, so it stays
    true whatever a rebuild did to the source, and it is deliberately
    conservative: a column that happens to hold only round numbers is reported
    as coarse, which can only make a tolerance harder to justify.
    """
    fields, grid = tables[output]
    if field not in fields:
        raise AcceptanceError(f"{output} has no {field} field")
    column = grid[:, :, fields.index(field)]
    for decimals in range(MAX_WRITTEN_DECIMALS + 1):
        if np.array_equal(np.round(column, decimals), column):
            return float(10.0 ** -decimals)
    return 0.0


def closure_resolution(tables: dict, contract: dict, count: int) -> dict:
    """What each closure residual can resolve, against the tolerance it is held to.

    A CONSERVATION CHECK CANNOT BE FINER THAN THE COLUMNS IT DIFFERENCES. The
    pool residual is `(pool[-1] - pool[0]) * factor + sum(flux)`: each pool
    endpoint carries up to half a written quantum, so their difference carries
    a whole one, and each of the summed flux terms carries half of its own. The
    water residual differences two summed columns against a driver
    precipitation that is full-precision binary and contributes nothing.

    The bound this returns is the residual an output written this coarsely
    reports for a model that conserves EXACTLY. A declared absolute floor at or
    below it is not a tolerance: every cell it fails, it fails for having been
    written down. `minimum_resolution_margin` is how far above its own
    quantisation a tolerance has to sit to be measuring the simulation.
    """
    closure = contract["closure"]
    margin = float(closure["minimum_resolution_margin"])
    resolution = {"minimum_resolution_margin": margin,
                  "form": "the residual an exactly conserving model reports, "
                          "given the written precision of the columns differenced"}
    unresolvable = []
    for element in ("carbon", "nitrogen"):
        rule = closure[element]
        factor = float(rule.get("pool_to_flux_units", 1.0))
        pool_q = written_quantum(tables, rule["pool_output"], rule["pool_field"])
        flux_q = written_quantum(tables, rule["flux_output"], rule["flux_field"])
        bound = pool_q * factor + (count - 1) * flux_q / 2.0
        floor = float(rule[absolute_floor_key(rule)])
        resolution[element] = {
            "pool_quantum": pool_q, "flux_quantum": flux_q,
            "pool_to_flux_units": factor,
            "resolution_bound": bound, "absolute_floor": floor,
            "margin": (floor / bound) if bound > 0 else None,
        }
        if bound > 0 and floor < margin * bound:
            unresolvable.append(
                f"{element} floor {floor:g} is {floor / bound:.2f}x its own "
                f"{bound:g} resolution bound, under the required {margin:g}x")
    water = closure["water"]
    aet_q = written_quantum(tables, water["aet_output"], water["aet_field"])
    runoff_q = written_quantum(tables, water["runoff_output"], water["runoff_field"])
    bound = count * (aet_q + runoff_q) / 2.0
    floor = float(water["absolute_floor_mm"])
    resolution["water"] = {
        "aet_quantum": aet_q, "runoff_quantum": runoff_q,
        "resolution_bound_mm": bound, "absolute_floor_mm": floor,
        "margin": (floor / bound) if bound > 0 else None,
    }
    if bound > 0 and floor < margin * bound:
        unresolvable.append(
            f"water floor {floor:g} mm is {floor / bound:.2f}x its own "
            f"{bound:g} mm resolution bound, under the required {margin:g}x")
    resolution["unresolvable"] = unresolvable
    return resolution


def absolute_floor_key(rule: dict) -> str:
    """The one `absolute_floor_*` key in a closure rule, named exactly once."""
    keys = [key for key in rule if key.startswith("absolute_floor_")]
    if len(keys) != 1:
        raise AcceptanceError(f"closure rule needs exactly one absolute floor, has {keys}")
    return keys[0]


def named_failures(cells: list[tuple[float, float]], residuals: np.ndarray,
                   limits: np.ndarray, failed: np.ndarray) -> list[dict]:
    """The failing cells a refusal names, largest residual first."""
    index = np.flatnonzero(failed)
    index = index[np.argsort(-np.abs(residuals[index]))][:NAMED_FAILURES]
    return [{"lon": float(cells[cell][0]), "lat": float(cells[cell][1]),
             "residual": float(residuals[cell]), "limit": float(limits[cell])}
            for cell in index]


def closure_report(tables: dict, cells: list[tuple[float, float]],
                   years: list[int], manifest: dict, contract: dict) -> dict:
    closure = contract["closure"]
    cycles = int(closure["window_complete_forcing_cycles"])
    cycle_years = int(manifest.get("forcing", {}).get("cycle_years", 1))
    count = cycles * cycle_years
    if len(years) < count:
        raise AcceptanceError(f"closure needs {count} end years")
    selected = range(len(years) - count, len(years))
    driver_record = manifest.get("inputs", {}).get("driver", {})
    driver = Path(driver_record.get("path", ""))
    if not driver.is_file() or sha256(driver) != driver_record.get("sha256"):
        raise AcceptanceError("manifest-pinned driver is absent or changed")
    precip, precip_identity = driver_precipitation(driver)
    if set(precip) != set(cells):
        raise AcceptanceError(
            f"driver/output support differs: {len(precip)} versus {len(cells)} cells")

    # THE INSTRUMENT BEFORE THE MEASUREMENT. Every residual below is a
    # difference of written columns, so a tolerance under the written precision
    # discriminates on rounding rather than on conservation. This refuses such a
    # tolerance by name instead of reporting its rounding as a defect.
    reports = {"resolution": closure_resolution(tables, contract, count)}
    if reports["resolution"]["unresolvable"]:
        raise ClosureError(
            "closure tolerance is finer than the output it reads: "
            + "; ".join(reports["resolution"]["unresolvable"]), reports)

    refusals = []
    for element in ("carbon", "nitrogen"):
        rule = closure[element]
        factor = float(rule.get("pool_to_flux_units", 1.0))
        floor = float(rule[absolute_floor_key(rule)])
        residuals = np.empty(len(cells))
        throughputs = np.empty(len(cells))
        for cell, _ in enumerate(cells):
            pool = field_series(tables, rule["pool_output"], rule["pool_field"],
                                cell, selected)
            flux = field_series(tables, rule["flux_output"], rule["flux_field"],
                                cell, selected)
            residuals[cell] = (pool[-1] - pool[0]) * factor + float(flux[1:].sum())
            throughputs[cell] = float(np.abs(flux[1:]).sum())
        limits = np.maximum(floor, float(rule["relative_throughput_limit"]) * throughputs)
        failed = np.abs(residuals) > limits
        worst = int(np.argmax(np.abs(residuals)))
        reports[element] = {
            "maximum_absolute_residual": float(np.abs(residuals[worst])),
            # The limit that applied to the WORST cell, which is the one the
            # maximum residual is to be read against. A maximum over every
            # cell's limit belongs to whichever cell had the largest
            # throughput and answers a question nobody asked.
            "limit_at_maximum_residual": float(limits[worst]),
            "absolute_floor": floor,
            "binding_limit": "absolute floor" if float(limits.max()) <= floor
                             else "mixed floor and relative throughput",
            "failed_cells": int(failed.sum()), "cells": len(cells),
            "named_failures": named_failures(cells, residuals, limits, failed),
        }
        if failed.any():
            refusals.append(f"{element} closure fails in {int(failed.sum())} cells")

    water = closure["water"]
    residuals = np.empty(len(cells))
    limits = np.empty(len(cells))
    for cell, coordinate in enumerate(cells):
        aet = field_series(tables, water["aet_output"], water["aet_field"],
                           cell, selected)
        runoff = field_series(tables, water["runoff_output"], water["runoff_field"],
                              cell, selected)
        p_total = cycles * float(precip[coordinate].sum())
        residuals[cell] = p_total - float(aet.sum()) - float(runoff.sum())
        limits[cell] = max(float(water["absolute_floor_mm"]),
                           float(water["relative_throughput_limit"]) * p_total)
    failed = np.abs(residuals) > limits
    worst = int(np.argmax(np.abs(residuals)))
    reports["water"] = {
        "form": "precipitation - AET - runoff = implied end-minus-start storage",
        "interpretation": water["interpretation"],
        "maximum_absolute_residual_mm": float(np.abs(residuals[worst])),
        "limit_at_maximum_residual_mm": float(limits[worst]),
        "failed_cells": int(failed.sum()), "cells": len(cells),
        "named_failures": named_failures(cells, residuals, limits, failed),
        "driver": precip_identity,
    }
    if failed.any():
        refusals.append(f"water closure fails in {int(failed.sum())} cells")

    # EVERY ELEMENT IS EVALUATED BEFORE ANY OF THEM REFUSES. Raising on carbon
    # left nitrogen and water unmeasured, so a run was fixed and re-run once per
    # element rather than once.
    if refusals:
        raise ClosureError("; ".join(refusals), reports)
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
        merged = read_table(run_dir / output)
        cells, years = coverage(merged, nyear)
        check_physical(output, merged.fields, merged.values, contract)
        pieces = []
        rank_cells = []
        for rank in range(1, ranks + 1):
            piece = read_table(run_dir / f"run{rank}" / output)
            if piece.fields != merged.fields:
                raise AcceptanceError(f"rank {rank} {output} has a different header")
            piece_cells, piece_years = coverage(piece, nyear)
            if piece_years != years:
                raise AcceptanceError(f"rank {rank} {output} has different years")
            pieces.append(piece)
            rank_cells.append(set(piece_cells))
        for left in range(ranks):
            for right in range(left + 1, ranks):
                if rank_cells[left] & rank_cells[right]:
                    raise AcceptanceError(f"{output} repeats cells across ranks")
        verify_merge(output, merged, pieces)
        counts = [len(value) for value in rank_cells]
        if max(counts) - min(counts) > 1:
            raise AcceptanceError(f"{output} rank cell counts are imbalanced: {counts}")
        if canonical_support is None:
            canonical_support, canonical_years = cells, years
        elif cells != canonical_support or years != canonical_years:
            raise AcceptanceError(f"{output} has different cell/year support")
        tables[output] = (merged.fields, merged.grid(len(cells), len(years)))
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
            # Both spans, because they are two different ones: `reported` is
            # what a consumer read the value over and `window` is the per-cell
            # half's alone.
            "reported": reduced.report.get("reported"),
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
    # A closure refusal knows which cells failed and by how much. Without this
    # the report carried a count, and a count is not something a reader can act
    # on: it says four cells are wrong and leaves finding them to whoever reads
    # the gigabytes again.
    if isinstance(error, ClosureError):
        report["closure"] = error.closure
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


# What the model writes each closure column to, from the ColumnDescriptors in
# `vendor/lpj-guess/modules/commonoutput.cpp`. The fixture writes at these
# precisions so the resolution check has a real quantum to measure; the check
# itself reads the quantum off the table and never off this map.
CLOSURE_DECIMALS = {
    ("cpool.out", "Total"): 6, ("cflux.out", "NEE"): 5,
    ("npool.out", "Total"): 7, ("nflux.out", "NEE"): 5,
    ("aaet.out", "Total"): 4, ("tot_runoff.out", "Total"): 4,
}


def _assessed_columns(output: str) -> list[str]:
    """Every column the equilibrium contract NAMES for one output table.

    The contract selects a quantity either by naming its columns or by naming
    the ones to exclude, and it refuses a table that carries neither. A fixture
    with a hand-written column list drifts the moment a quantity is added to
    the contract, and then the self-test exercises that drift refusal instead
    of the checks it was written for. So the fixture's columns come from the
    contract, and adding a row there needs no edit here.
    """
    policy = yaml.safe_load(EQUILIBRIUM_CONFIG.read_text(encoding="utf-8"))
    named = []
    for quantity in policy.get("assessed", {}).get("quantities", []):
        if quantity.get("table") != output:
            continue
        for name in (quantity.get("columns") or []) + (quantity.get("columns_excluding") or []):
            if name not in named:
                named.append(name)
    return named


def _table_text(output: str, cells: list[tuple[float, float]], years: range,
                coarsen: dict | None = None) -> str:
    """One fixture output table, written the way the model writes it.

    Each closure column varies in its LAST written digit from year to year. A
    column holding one constant reads back as coarser than it was written --
    every value in it is round at any precision -- so a constant fixture would
    exercise the resolution check's degenerate case instead of its ordinary one.
    The variation is a single quantum, which is negligible against every
    tolerance in the contract and keeps the fixture a passing run.
    """
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
    # Any column the equilibrium contract names that the fixture does not
    # already carry, so the contract and the fixture cannot drift apart.
    for name in _assessed_columns(output):
        if name not in fields:
            fields, values = fields + [name], values + [0.1]
    decimals = [(coarsen or {}).get(field, CLOSURE_DECIMALS.get((output, field)))
                for field in fields]
    lines = ["Lon Lat Year " + " ".join(fields)]
    for lon, lat in cells:
        for year in years:
            written = []
            for value, places in zip(values, decimals):
                if places is None:
                    written.append(str(value))
                else:
                    written.append(f"{value + (year % 10) * 10.0 ** -places:.{places}f}")
            lines.append(f"{lon} {lat} {year} " + " ".join(written))
    return "\n".join(lines) + "\n"


def _parsers_agree(root: Path) -> bool:
    """Do the column-wise and row-at-a-time parsers read one table the same?

    The whole reason the gate can read a table column-wise is that the two
    parsers are the same parser, so their agreement is a fixture rather than an
    assumption. The bed is deliberately awkward for the column-wise path:
    coordinates at both poles and either edge of the antimeridian, negative and
    zero values, and exponents.
    """
    path = root / "parser_agreement.out"
    path.write_text(
        "Lon Lat Year A B\n"
        "-179.95 -89.75 4 1.5 -2.25\n"
        "-179.95 -89.75 3 0.0 1.0e-8\n"
        "  0.05  89.75 4 -0.5 3.0\n"
        "  0.05  89.75 3 2.0 -1.0e+10\n"
        " 179.95   0.00 4 1.0 2.0\n"
        " 179.95   0.00 3 3.0 4.0\n", encoding="utf-8")
    fast = lpj_table.read(path)
    slow = read_table_rows(path)
    return (fast.fields == slow.fields
            and np.array_equal(fast.lon, slow.lon)
            and np.array_equal(fast.lat, slow.lat)
            and np.array_equal(fast.year, slow.year)
            and np.array_equal(fast.key, slow.key)
            and np.array_equal(fast.values, slow.values))


def selftest() -> dict:
    fixtures = []
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        run = root / "run"
        run.mkdir()
        cells = [(0.0, 10.0), (20.0, -10.0)]
        # The equilibrium contract measures its own trending-cell null on the
        # record BEFORE the acceptance window, so a fixture has to retain one:
        # at a per-field false-refusal rate of 0.05 that is 19 null windows
        # beside the acceptance window, and 20 windows of ten one-year cycles is
        # 200 years plus the window. A fixture shorter than that exercises the
        # insufficient-record refusal instead of the checks it was written for.
        years = range(210)
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
            "physical": {"nyear": len(years), "ranks": 2}, "ranks": 2,
            "forcing": {"cells": 2, "cycle_years": 1},
            "inputs": {"driver": {"path": str(driver), "sha256": sha256(driver)}},
        }
        (run / "run_manifest.json").write_text(json.dumps(manifest) + "\n")
        assess(run, write=False)
        fixtures.append({"fixture": "complete run", "pass": True})
        fixtures.append({"fixture": "column and row parsers agree",
                         "pass": _parsers_agree(root)})

        def coarsen_npool(bed: Path) -> None:
            """Write npool.out at the precision the release shipped.

            Four decimals in kgN/m2 against a `pool_to_flux_units` of 1e4 is one
            least-significant digit per 1.0 kgN/ha, so the 2.0 kgN/ha floor sits
            at twice its own quantisation and fails cells for how they were
            written down. The gate has to refuse the TOLERANCE here, not the
            cells, and it has to say so.
            """
            for rank, cell in enumerate(cells, 1):
                (bed / f"run{rank}" / "npool.out").write_text(
                    _table_text("npool.out", [cell], years, coarsen={"Total": 4}),
                    encoding="utf-8")
            parts = [(bed / f"run{rank}" / "npool.out").read_text(
                encoding="utf-8").splitlines() for rank in (1, 2)]
            (bed / "npool.out").write_text(
                "\n".join(parts[0] + parts[1][1:]) + "\n", encoding="utf-8")

        mutations = {
            "missing rank output": (
                lambda bed: (bed / "run2" / "lai.out").unlink(), None),
            "non-finite value": (
                lambda bed: (bed / "anpp.out").write_text(
                    (bed / "anpp.out").read_text().replace("1.0", "nan", 1)), None),
            "negative stock": (
                lambda bed: (bed / "cpool.out").write_text(
                    (bed / "cpool.out").read_text().replace("10.0", "-1.0", 1)), None),
            "carbon leak": (
                lambda bed: (bed / "cflux.out").write_text(
                    (bed / "cflux.out").read_text().replace("0.0", "1.0")), None),
            "water leak": (
                lambda bed: (bed / "aaet.out").write_text(
                    (bed / "aaet.out").read_text().replace("600.0", "500.0")), None),
            "closure tolerance under the written precision": (
                coarsen_npool, "finer than the output it reads"),
        }
        for name, (mutate, expected) in mutations.items():
            bed = root / name.replace(" ", "_")
            shutil.copytree(run, bed)
            mutate(bed)
            try:
                assess(bed, write=False)
            except AcceptanceError as exc:
                # A mutation that names its refusal has to earn THAT refusal.
                # Any AcceptanceError would otherwise pass a fixture whose whole
                # point is which check fired.
                fixtures.append({"fixture": name,
                                 "pass": expected is None or expected in str(exc)})
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
