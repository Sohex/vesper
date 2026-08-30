"""One equilibrium-window reducer for every LPJ-GUESS output consumer."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import struct
from typing import Iterable

import numpy as np
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
POLICY_PATH = PROJECT_ROOT / "biosphere/config/equilibrium_window.yaml"
DRIVER_MAGIC = b"VESPDRV8"


class EquilibriumWindowError(ValueError):
    """An LPJ table cannot support the declared equilibrium statistic."""


def require_lpj_acceptance(path: Path) -> dict:
    """Require BIO-14's PASS artifact and prove it still covers this input."""
    path = Path(path)
    run_dir = path if path.is_dir() else path.parent
    report_path = run_dir / "acceptance.json"
    manifest_path = run_dir / "run_manifest.json"
    if not report_path.is_file():
        raise EquilibriumWindowError(
            f"{path} has no sibling BIO-14 acceptance.json")
    try:
        report = json.loads(report_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise EquilibriumWindowError(f"cannot read {report_path}: {exc}") from exc
    if report.get("verdict") != "PASS":
        raise EquilibriumWindowError(
            f"{report_path} refuses this run: {report.get('refusal', 'no reason recorded')}")
    if not manifest_path.is_file() or sha256(manifest_path) != report.get("manifest_sha256"):
        raise EquilibriumWindowError(
            f"{manifest_path} is absent or changed since BIO-14 assessment")
    if path.is_file():
        expected = report.get("output_sha256", {}).get(path.name)
        if expected is None or sha256(path) != expected:
            raise EquilibriumWindowError(
                f"{path} is absent from or changed since BIO-14 assessment")
    return report


@dataclass
class ReducedTable:
    names: list[str]
    values: dict[tuple[float, float], np.ndarray]
    temporal_std: dict[tuple[float, float], np.ndarray]
    report: dict


def read_policy(path: Path = POLICY_PATH) -> dict:
    policy = yaml.safe_load(path.read_text())
    if policy.get("contract_version") != "vesper-lpj-equilibrium-window/1":
        raise EquilibriumWindowError("unsupported equilibrium-window contract")
    cycles = policy.get("complete_forcing_cycles")
    if not isinstance(cycles, int) or cycles < 3:
        raise EquilibriumWindowError("complete_forcing_cycles must be at least 3")
    trend = policy.get("trend", {})
    for key in ("relative_end_to_end_limit", "slope_standard_errors",
                "maximum_trending_cell_fraction", "absolute_scale_floor"):
        if not isinstance(trend.get(key), (int, float)) or trend[key] <= 0:
            raise EquilibriumWindowError(f"trend.{key} must be positive")
    return policy


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def driver_header(path: Path) -> dict:
    with path.open("rb") as handle:
        magic = handle.read(8)
        raw = handle.read(24)
    if magic != DRIVER_MAGIC or len(raw) != 24:
        raise EquilibriumWindowError(
            f"{path} is not a complete {DRIVER_MAGIC.decode()} driver")
    ncells, intervals, year_days, years, subdaily, padding = struct.unpack(
        "<6i", raw)
    if ncells < 1 or intervals < 1 or year_days < 1 or years < 1:
        raise EquilibriumWindowError(f"{path} has an invalid driver header")
    return {"format": magic.decode(), "cells": ncells,
            "intervals_per_year": intervals, "year_length_days": year_days,
            "cycle_years": years, "subdaily_samples": subdaily,
            "header_padding": padding}


def _manifest_for(path: Path) -> tuple[Path, dict]:
    manifest_path = path.parent / "run_manifest.json"
    if not manifest_path.is_file():
        raise EquilibriumWindowError(
            f"{path} has no sibling run_manifest.json; forcing-cycle length "
            "must be provenance, never inferred from the output")
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise EquilibriumWindowError(f"cannot read {manifest_path}: {exc}") from exc
    return manifest_path, manifest


def forcing_cycle_years(path: Path, manifest: dict) -> tuple[int, dict]:
    forcing = manifest.get("forcing", {})
    cycle = forcing.get("cycle_years")
    if isinstance(cycle, int) and cycle >= 1:
        return cycle, {"source": "run_manifest.forcing", **forcing}

    driver = manifest.get("inputs", {}).get("driver", {})
    driver_path = driver.get("path")
    expected_hash = driver.get("sha256")
    if not driver_path:
        raise EquilibriumWindowError(
            f"{path}'s manifest carries neither forcing.cycle_years nor a driver path")
    driver_path = Path(driver_path)
    if not driver_path.is_file():
        raise EquilibriumWindowError(
            f"{path}'s forcing-cycle source {driver_path} is absent; regenerate "
            "the run manifest with an embedded forcing header")
    actual_hash = sha256(driver_path)
    if expected_hash and actual_hash != expected_hash:
        raise EquilibriumWindowError(
            f"{driver_path} changed after the run; its hash cannot establish "
            f"{path}'s forcing cycle")
    header = driver_header(driver_path)
    return header["cycle_years"], {
        "source": "manifest-pinned driver header", "driver": str(driver_path),
        "driver_sha256": actual_hash, **header}


def _read_rows(path: Path) -> tuple[list[str], list[tuple[float, float]],
                                    list[int], np.ndarray]:
    lines = path.read_text().splitlines()
    if not lines:
        raise EquilibriumWindowError(f"{path} is empty")
    header = lines[0].split()
    if len(header) < 4 or header[:3] != ["Lon", "Lat", "Year"]:
        raise EquilibriumWindowError(
            f"{path}'s first columns must be exactly Lon Lat Year")
    names = header[3:]
    rows: dict[tuple[tuple[float, float], int], np.ndarray] = {}
    for line_number, line in enumerate(lines[1:], 2):
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) != len(header):
            raise EquilibriumWindowError(
                f"{path}:{line_number} has {len(parts)} fields, expected {len(header)}")
        try:
            lon, lat = round(float(parts[0]), 2), round(float(parts[1]), 2)
            year_value = float(parts[2])
            year = int(year_value)
            values = np.asarray([float(value) for value in parts[3:]], dtype=float)
        except ValueError as exc:
            raise EquilibriumWindowError(
                f"{path}:{line_number} contains a nonnumeric field") from exc
        if year_value != year:
            raise EquilibriumWindowError(f"{path}:{line_number} has nonintegral year")
        if not np.isfinite(values).all():
            raise EquilibriumWindowError(f"{path}:{line_number} contains NaN or infinity")
        identity = ((lon, lat), year)
        if identity in rows:
            raise EquilibriumWindowError(
                f"{path}:{line_number} duplicates cell {(lon, lat)} year {year}")
        rows[identity] = values
    if not rows:
        raise EquilibriumWindowError(f"{path} contains no data rows")

    cells = sorted({cell for cell, _ in rows})
    maxima: dict[tuple[float, float], int] = {}
    for cell, year in rows:
        maxima[cell] = max(year, maxima.get(cell, year))
    if len(set(maxima.values())) != 1:
        detail = sorted(set(maxima.values()))
        raise EquilibriumWindowError(
            f"{path} ended at different years across cells: {detail}")
    latest = next(iter(maxima.values()))
    years = sorted({year for _, year in rows})
    cube = np.full((len(years), len(cells), len(names)), np.nan)
    year_index = {year: index for index, year in enumerate(years)}
    cell_index = {cell: index for index, cell in enumerate(cells)}
    for (cell, year), values in rows.items():
        cube[year_index[year], cell_index[cell]] = values
    return names, cells, years, cube


def _trend(cycle_means: np.ndarray, policy: dict) -> tuple[list[dict], bool]:
    # cycle_means: cycle, cell, field
    ncycle, ncell, nfield = cycle_means.shape
    x = np.arange(ncycle, dtype=float)
    x -= x.mean()
    denominator = float(np.sum(x * x))
    slope = np.einsum("t,tcf->cf", x, cycle_means) / denominator
    intercept = cycle_means.mean(axis=0)
    residual = cycle_means - (intercept[None] + x[:, None, None] * slope[None])
    if ncycle > 2:
        slope_se = np.sqrt(np.sum(residual * residual, axis=0)
                           / (ncycle - 2) / denominator)
    else:  # read_policy requires at least three, retained as a hard backstop.
        slope_se = np.full_like(slope, np.inf)
    drift = slope * (ncycle - 1)
    floor = float(policy["trend"]["absolute_scale_floor"])
    scale = np.maximum(np.abs(intercept), floor)
    relative = np.abs(drift) / scale
    sigma = float(policy["trend"]["slope_standard_errors"])
    with np.errstate(divide="ignore", invalid="ignore"):
        significance = np.where(slope_se > 0, np.abs(slope) / slope_se,
                                np.where(np.abs(slope) > 0, np.inf, 0.0))
    cell_trending = ((relative > policy["trend"]["relative_end_to_end_limit"])
                     & (significance > sigma))

    spatial = cycle_means.mean(axis=1)
    spatial_slope = np.einsum("t,tf->f", x, spatial) / denominator
    spatial_mean = spatial.mean(axis=0)
    spatial_residual = spatial - (spatial_mean[None]
                                  + x[:, None] * spatial_slope[None])
    spatial_se = np.sqrt(np.sum(spatial_residual * spatial_residual, axis=0)
                         / (ncycle - 2) / denominator)
    spatial_drift = spatial_slope * (ncycle - 1)
    spatial_relative = np.abs(spatial_drift) / np.maximum(np.abs(spatial_mean), floor)
    with np.errstate(divide="ignore", invalid="ignore"):
        spatial_significance = np.where(
            spatial_se > 0, np.abs(spatial_slope) / spatial_se,
            np.where(np.abs(spatial_slope) > 0, np.inf, 0.0))
    global_trending = (
        (spatial_relative > policy["trend"]["relative_end_to_end_limit"])
        & (spatial_significance > sigma))
    fractions = cell_trending.sum(axis=0) / max(ncell, 1)
    fraction_trending = (
        fractions > policy["trend"]["maximum_trending_cell_fraction"])
    rejected = global_trending | fraction_trending
    diagnostics = [{
        "spatial_mean_relative_end_to_end_change": float(spatial_relative[i]),
        "spatial_mean_slope_standard_errors": float(spatial_significance[i]),
        "trending_cell_fraction": float(fractions[i]),
        "global_trending": bool(global_trending[i]),
        "cell_fraction_trending": bool(fraction_trending[i]),
        "rejected": bool(rejected[i]),
    } for i in range(nfield)]
    return diagnostics, bool(rejected.any())


def _summary(values: np.ndarray, names: list[str]) -> dict:
    # values: cell, field
    return {name: {
        "mean": float(np.mean(values[:, i])),
        "median": float(np.median(values[:, i])),
        "p95": float(np.quantile(values[:, i], 0.95)),
        "maximum": float(np.max(values[:, i])),
    } for i, name in enumerate(names)}


def _physical_identity(manifest: dict) -> dict:
    physical = dict(manifest.get("physical", {}))
    for key in ("root_seed", "npatch", "ranks", "label"):
        physical.pop(key, None)
    inputs = manifest.get("inputs", {})
    return {
        "source_build": manifest.get("source_build"),
        "physical_except_sampling": physical,
        "inputs": {name: value.get("sha256") for name, value in inputs.items()
                   if name != "binary" and isinstance(value, dict)},
    }


def reduce_table(path: Path, peers: Iterable[Path] = (), *,
                 policy_path: Path = POLICY_PATH) -> ReducedTable:
    path = Path(path)
    policy = read_policy(policy_path)
    manifest_path, manifest = _manifest_for(path)
    cycle_years, forcing_source = forcing_cycle_years(path, manifest)
    names, cells, years, cube = _read_rows(path)
    window_years = policy["complete_forcing_cycles"] * cycle_years
    latest = max(years)
    selected_years = list(range(latest - window_years + 1, latest + 1))
    if len(years) < window_years or any(year not in years for year in selected_years):
        raise EquilibriumWindowError(
            f"{path} needs {window_years} consecutive end years "
            f"({policy['complete_forcing_cycles']} complete {cycle_years}-year "
            f"forcing cycles), ending at {latest}")
    selected_indices = [years.index(year) for year in selected_years]
    window = cube[selected_indices]
    if not np.isfinite(window).all():
        missing = int(np.size(window) - np.isfinite(window).sum())
        raise EquilibriumWindowError(
            f"{path}'s declared end window has {missing} missing cell-year rows")
    mean = window.mean(axis=0)
    temporal_std = window.std(axis=0, ddof=1)
    cycle_means = window.reshape(
        policy["complete_forcing_cycles"], cycle_years, len(cells), len(names)
    ).mean(axis=1)
    trends, rejected = _trend(cycle_means, policy)
    for name, diagnostic in zip(names, trends):
        diagnostic["field"] = name
    if rejected:
        failed = [item["field"] for item in trends if item["rejected"]]
        raise EquilibriumWindowError(
            f"{path}'s end window is still trending in {', '.join(failed)}")

    peer_records = []
    ensembles = [(path, manifest, mean)]
    identity = _physical_identity(manifest)
    for peer_path in peers:
        peer = reduce_table(Path(peer_path), policy_path=policy_path)
        _, peer_manifest = _manifest_for(Path(peer_path))
        if peer.names != names or sorted(peer.values) != cells:
            raise EquilibriumWindowError(f"{peer_path} has a different table support")
        if _physical_identity(peer_manifest) != identity:
            raise EquilibriumWindowError(
                f"{peer_path} differs in physical inputs, not only seed/patch sampling")
        peer_mean = np.stack([peer.values[cell] for cell in cells])
        ensembles.append((Path(peer_path), peer_manifest, peer_mean))
        peer_records.append(peer.report["identity"])

    sampling = []
    for member_path, member_manifest, member_mean in ensembles:
        physical = member_manifest.get("physical", {})
        stochastic = member_manifest.get("stochastic_randomness", {})
        sampling.append({"path": str(member_path),
                         "root_seed": stochastic.get(
                             "root_seed", physical.get("root_seed")),
                         "npatch": physical.get("npatch"), "mean": member_mean})
    uncertainty: dict = {
        "temporal": {
            "statistic": policy["uncertainty"]["temporal_statistic"],
            "by_field_across_cells": _summary(temporal_std, names),
        },
        "seed": {"status": "not_measured",
                 "reason": "fewer than two root seeds at one patch count"},
        "patch_count": {"status": "not_measured",
                        "reason": "fewer than two patch counts"},
        "ensemble_members": [{k: v for k, v in item.items() if k != "mean"}
                             for item in sampling],
    }
    by_patch: dict[int, list[dict]] = {}
    for item in sampling:
        if isinstance(item["npatch"], int):
            by_patch.setdefault(item["npatch"], []).append(item)
    seed_groups = [group for group in by_patch.values()
                   if len({item["root_seed"] for item in group
                           if item["root_seed"] is not None}) >= 2]
    if seed_groups:
        spread = np.mean([
            np.stack([item["mean"] for item in group]).std(axis=0, ddof=1)
            for group in seed_groups], axis=0)
        uncertainty["seed"] = {
            "status": "measured",
            "statistic": policy["uncertainty"]["seed_statistic"],
            "by_field_across_cells": _summary(spread, names),
            "patch_counts": sorted(group[0]["npatch"] for group in seed_groups),
        }
    if len(by_patch) >= 2:
        patch_means = np.stack([
            np.stack([item["mean"] for item in group]).mean(axis=0)
            for _, group in sorted(by_patch.items())])
        spread = patch_means.max(axis=0) - patch_means.min(axis=0)
        uncertainty["patch_count"] = {
            "status": "measured",
            "statistic": policy["uncertainty"]["patch_statistic"],
            "by_field_across_cells": _summary(spread, names),
            "patch_counts": sorted(by_patch),
        }

    report = {
        "contract_version": policy["contract_version"],
        "policy": str(policy_path.relative_to(PROJECT_ROOT)),
        "policy_sha256": sha256(policy_path),
        "identity": {"table": str(path), "run_manifest": str(manifest_path),
                     "run_id": manifest.get("run_id"),
                     "root_seed": manifest.get("stochastic_randomness", {}).get(
                         "root_seed", manifest.get("physical", {}).get("root_seed")),
                     "npatch": manifest.get("physical", {}).get("npatch")},
        "forcing": forcing_source,
        "window": {"first_year": selected_years[0], "last_year": selected_years[-1],
                   "annual_values": window_years,
                   "complete_forcing_cycles": policy["complete_forcing_cycles"],
                   "forcing_cycle_years": cycle_years},
        "trend": {"rule": policy["trend"], "fields": trends, "verdict": "PASS"},
        "uncertainty": uncertainty,
        "peers": peer_records,
    }
    return ReducedTable(
        names=names,
        values={cell: mean[i] for i, cell in enumerate(cells)},
        temporal_std={cell: temporal_std[i] for i, cell in enumerate(cells)},
        report=report)
