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

from autocorrelation import (RELIABLE_SPAN_MULTIPLE, integrated_time,
                             mean_standard_error)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
POLICY_PATH = PROJECT_ROOT / "biosphere/config/equilibrium_window.yaml"
MEMORY_ESTIMATOR = "lib/autocorrelation.py:integrated_time"
ACCEPTANCE_PATH = PROJECT_ROOT / "biosphere/config/lpj_acceptance.yaml"
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
    if policy.get("contract_version") != "vesper-lpj-equilibrium-window/3":
        raise EquilibriumWindowError("unsupported equilibrium-window contract")
    cycles = policy.get("complete_forcing_cycles")
    if not isinstance(cycles, int) or cycles < 3:
        raise EquilibriumWindowError("complete_forcing_cycles must be at least 3")
    trend = policy.get("trend", {})
    for key in ("relative_end_to_end_limit", "slope_standard_errors",
                "absolute_scale_floor"):
        if not isinstance(trend.get(key), (int, float)) or trend[key] <= 0:
            raise EquilibriumWindowError(f"trend.{key} must be positive")
    cell = trend.get("cell_fraction", {})
    if cell.get("null_statistic") != "maximum_over_detrended_windows_of_this_run":
        raise EquilibriumWindowError(
            "trend.cell_fraction.null_statistic names the only null this "
            "contract knows how to build")
    rate = cell.get("per_field_false_refusal_rate")
    if not isinstance(rate, (int, float)) or not 0 < rate < 1:
        raise EquilibriumWindowError(
            "trend.cell_fraction.per_field_false_refusal_rate must lie in (0, 1)")
    if trend.get("memory_estimator") != MEMORY_ESTIMATOR:
        raise EquilibriumWindowError(
            f"trend.memory_estimator must name {MEMORY_ESTIMATOR}")
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


def _cell_fraction(cycle_means: np.ndarray, policy: dict) -> np.ndarray:
    """Per field, the share of the cells it OCCUPIES whose own series is trending.

    The denominator is the occupied cells and not every cell, so the statistic
    separates how widely a field drifts from how widely it is present. Over the
    whole grid a plant functional type on a sixth of the cells could not reach a
    quarter-of-all-cells limit however hard it drifted, which made the test inert
    for nine of this world's thirteen types.
    """
    ncycle = cycle_means.shape[0]
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
    floor = float(policy["trend"]["absolute_scale_floor"])
    relative = np.abs(slope * (ncycle - 1)) / np.maximum(np.abs(intercept), floor)
    with np.errstate(divide="ignore", invalid="ignore"):
        significance = np.where(slope_se > 0, np.abs(slope) / slope_se,
                                np.where(np.abs(slope) > 0, np.inf, 0.0))
    trending = ((relative > policy["trend"]["relative_end_to_end_limit"])
                & (significance > float(policy["trend"]["slope_standard_errors"])))
    occupied = np.abs(intercept) > floor
    counts = occupied.sum(axis=0)
    return np.where(counts > 0, (trending & occupied).sum(axis=0)
                    / np.maximum(counts, 1), 0.0)


def _trend(cycle_means: np.ndarray, policy: dict,
           limits: np.ndarray) -> tuple[list[dict], bool]:
    # cycle_means: cycle, cell, field
    ncycle, ncell, nfield = cycle_means.shape
    x = np.arange(ncycle, dtype=float)
    x -= x.mean()
    denominator = float(np.sum(x * x))
    floor = float(policy["trend"]["absolute_scale_floor"])
    sigma = float(policy["trend"]["slope_standard_errors"])

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
    fractions = _cell_fraction(cycle_means, policy)
    fraction_trending = fractions > limits
    rejected = global_trending | fraction_trending
    diagnostics = [{
        "spatial_mean_relative_end_to_end_change": float(spatial_relative[i]),
        "spatial_mean_slope_standard_errors": float(spatial_significance[i]),
        "trending_cell_fraction": float(fractions[i]),
        "trending_cell_fraction_limit": float(limits[i]),
        "global_trending": bool(global_trending[i]),
        "cell_fraction_trending": bool(fraction_trending[i]),
        "rejected": bool(rejected[i]),
    } for i in range(nfield)]
    return diagnostics, bool(rejected.any())


def relaxation_time(series: np.ndarray, tau_memory: float,
                    blocks: int = 4) -> dict:
    """The e-folding time of an approach, measured WITHOUT its asymptote.

    WHY NOT A CURVE FIT. Fitting `a + b * exp(-t / tau)` needs the record to
    contain the turn-over: the asymptote `a` is a free parameter, and on a record
    shorter than the approach it lands outside the data and the fit says nothing.
    On this model's 1000-cycle record that happened for 34 of 64 assessed fields.

    WHAT THIS DOES INSTEAD. For that same exponential the DIFFERENCE between
    consecutive equal blocks decays by `exp(-Q / tau)`, and the asymptote cancels
    out of the ratio. So the record is cut into `blocks` equal parts, and the ratio
    of successive differences of their means estimates the decay RATE rather than
    the distance still to travel. A record shorter than tau can carry that.

    IT FAILS SOFTLY, so every condition it has to clear is declared rather than
    discovered: the successive differences must share a sign, because opposite
    signs are not a monotone approach; their ratio must be a contraction in (0, 1),
    because a series that is not shrinking has no e-folding time; each difference
    must exceed twice its own standard error, and that error is the
    MEMORY-CORRECTED one, because these block means carry the memory this module
    exists to respect; and the second ratio must agree with the first inside a
    factor of two, because one exponential has one rate. A field failing any of
    them gets NO relaxation time and is not given a default.
    """
    n = int(series.size)
    q = n // blocks
    if q < 3:
        return {"admissible": False, "reason": "fewer than three cycles per block"}
    parts = [series[i * q:(i + 1) * q] for i in range(blocks)]
    means = [float(part.mean()) for part in parts]
    errors = [float(mean_standard_error(part, tau_memory)) for part in parts]
    steps = [means[i + 1] - means[i] for i in range(blocks - 1)]
    step_errors = [float(np.hypot(errors[i], errors[i + 1]))
                   for i in range(blocks - 1)]
    base = {"block_cycles": q, "steps": steps, "step_standard_errors": step_errors}
    if not steps[0] or not steps[1] or (steps[0] > 0) != (steps[1] > 0):
        return {**base, "admissible": False,
                "reason": "successive block differences change sign, so the "
                          "series is not a monotone approach"}
    for index in (0, 1):
        if abs(steps[index]) <= 2.0 * step_errors[index]:
            return {**base, "admissible": False,
                    "reason": "a block difference is inside twice its own "
                              "memory-corrected standard error, so the approach "
                              "is smaller than the instrument reading it"}
    ratio = steps[1] / steps[0]
    if not 0.0 < ratio < 1.0:
        return {**base, "admissible": False, "ratio": ratio,
                "reason": f"a ratio of {ratio:.3f} is not a contraction"}
    tau = -q / float(np.log(ratio))
    second = None
    if (steps[2] and (steps[2] > 0) == (steps[1] > 0)
            and abs(steps[2]) > 2.0 * step_errors[2]):
        check = steps[2] / steps[1]
        if 0.0 < check < 1.0:
            second = -q / float(np.log(check))
            if not 0.5 <= second / tau <= 2.0:
                return {**base, "admissible": False, "ratio": ratio,
                        "tau_cycles": tau, "second_estimate_cycles": second,
                        "reason": f"the two ratios give {tau:.0f} and "
                                  f"{second:.0f} cycles, further than a factor "
                                  "of two apart, so this is not one exponential"}
    return {**base, "admissible": True, "ratio": ratio, "tau_cycles": tau,
            "second_estimate_cycles": second}


def timescale_report(run_dir: Path, *, policy_path: Path = POLICY_PATH) -> dict:
    """Both ecological timescales, per field, for every assessed table.

    Recorded on a run's acceptance artifact whatever its verdict, because a run
    that is refused for not having settled is exactly the run whose timescales say
    how long the next one has to be. `lib/run_lengths.py` reads this and states no
    number of its own, on the same terms as the climate relaxation bracket.
    """
    run_dir = Path(run_dir)
    policy = read_policy(policy_path)
    floor = float(policy["trend"]["absolute_scale_floor"])
    outputs = yaml.safe_load(ACCEPTANCE_PATH.read_text())["stability_outputs"]
    tables = {}
    for output in outputs:
        path = run_dir / output
        if not path.is_file():
            continue
        _, manifest = _manifest_for(path)
        cycle_years, _ = forcing_cycle_years(path, manifest)
        names, cells, years, cube = _read_rows(path)
        usable = (len(years) // cycle_years) * cycle_years
        record = cube[:usable].reshape(
            usable // cycle_years, cycle_years, len(cells), len(names)).mean(axis=1)
        span = record.shape[0]
        x = np.arange(span, dtype=float)
        fields = {}
        for index, name in enumerate(names):
            spatial = record[:, :, index].mean(axis=1)
            if (abs(float(spatial.mean())) <= floor
                    or not np.isfinite(spatial).all() or spatial.std() == 0):
                continue
            flat = spatial - np.polyval(np.polyfit(x, spatial, 1), x) + spatial.mean()
            memory = integrated_time(flat)
            fields[name] = {
                "memory": {"tau_cycles": float(memory["tau"]),
                           "effective_samples": float(memory["effective_sample_size"]),
                           "reliable": bool(memory["reliable"]),
                           "lag1": float(memory["lag1"])},
                "relaxation": relaxation_time(spatial, memory["tau"]),
            }
        tables[output] = {"record_cycles": int(span),
                          "forcing_cycle_years": int(cycle_years),
                          "fields": fields}
    return {"estimator": MEMORY_ESTIMATOR,
            "relaxation_estimator": "lib/lpj_output.py:relaxation_time",
            "reliable_span_multiple": float(RELIABLE_SPAN_MULTIPLE),
            # Carried so that `lib/run_lengths.py` can size a spin-up from this
            # artifact alone and state no number of its own. The residual a
            # spin-up has to decay to is the drift the contract that judges what
            # follows will absorb, which is the same reasoning
            # `SETTLING_RESIDUAL_K` uses on the climate side.
            "drift_tolerance": float(policy["trend"]["relative_end_to_end_limit"]),
            "contract_version": policy["contract_version"],
            "tables": tables}


def _memory_adequacy(record: np.ndarray, names: list[str],
                     policy: dict) -> tuple[list[dict], list[str]]:
    """Can this retained record establish the memory time of what it is judging?

    THE MEASUREMENT THAT FORCED THIS. At fixed forcing the integrated
    autocorrelation time of these spatial-mean series is 9 to 125 complete forcing
    cycles, against an acceptance window of 10. A trend fitted inside one memory
    time is not a trend: it is one smooth excursion of a process that has not had
    time to sample its own distribution, its residuals are small because the
    process is smooth on that scale, and the ordinary least-squares standard error
    it is judged against is correspondingly small. That is the whole of the 0.47 to
    0.86 per-cell flag rate this contract measures against a nominal 0.081.
    `biosphere/notes/equilibrium-trend-null.md` carries the measurement.

    WHAT THIS GUARD DOES, AND WHAT IT DOES NOT. It applies
    `lib/autocorrelation.py`'s declared span bar, fixed before it was applied to
    any series here, to BOTH spans that have to carry the memory time: the retained
    record, which is where tau is estimated, and the acceptance window, which is
    where the trend test below actually runs. A statistic computed over a span
    shorter than that span's own memory time is not supported by it, and the window
    is the span the slope is fitted on.

    It does not repair the trend test. That test is valid only where its window is
    long against tau; on this model today it is nowhere, so this guard fails closed
    and the test is unreachable. Its replacement is tracked rather than improvised.
    """
    floor = float(policy["trend"]["absolute_scale_floor"])
    window = int(policy["complete_forcing_cycles"])
    diagnostics, refused = [], []
    span = record.shape[0]
    x = np.arange(span, dtype=float)
    for index, name in enumerate(names):
        spatial = record[:, :, index].mean(axis=1)
        level = float(abs(spatial.mean()))
        if level <= floor or not np.isfinite(spatial).all() or spatial.std() == 0:
            diagnostics.append({
                "field": name, "assessed": False,
                "reason": "the spatial mean is zero or exactly constant, so it "
                          "carries no memory to establish"})
            continue
        flat = spatial - np.polyval(np.polyfit(x, spatial, 1), x) + spatial.mean()
        memory = integrated_time(flat)
        needed = float(RELIABLE_SPAN_MULTIPLE * memory["tau"])
        row = {"field": name, "assessed": True, "tau_cycles": float(memory["tau"]),
               "effective_samples": float(memory["effective_sample_size"]),
               "record_supports_tau": bool(memory["reliable"]),
               "window_supports_a_trend": bool(window >= needed),
               "cycles_required": needed}
        diagnostics.append(row)
        if not memory["reliable"]:
            refused.append(
                f"{name} (a {span}-cycle record cannot establish a memory time of "
                f"{memory['tau']:.1f} cycles, which needs {needed:.0f})")
        elif window < needed:
            refused.append(
                f"{name} (memory time {memory['tau']:.1f} cycles, so a trend over "
                f"the {window}-cycle window is inside one memory time; it needs "
                f"{needed:.0f})")
    return diagnostics, refused


def _cell_fraction_null(cube: np.ndarray, years: list[int], window_years: int,
                        cycle_years: int, policy: dict) -> tuple[np.ndarray, dict]:
    """Per field, the trending-cell fraction this run reaches with no trend left.

    The limit on the trending-cell fraction cannot be one declared number. Measured
    on a run at fixed forcing, the fraction a field reaches when nothing is drifting
    spans three orders of magnitude across the assessed fields, because it is set by
    that field's own internal variability and not by anything about equilibrium: the
    slow soil pools sit near a thousandth while the patch-driven grass and vegetation
    fields sit near four tenths. A single limit is therefore either unreachable for
    one group or permanently exceeded by the other.

    So the limit is measured rather than declared, and from the run being judged. The
    record before the acceptance window is detrended cell by cell, which removes any
    approach to equilibrium the run does carry, and is cut into windows of the
    contract's own length. A field's limit is the largest trending-cell fraction it
    reaches over those windows. The acceptance window itself is judged undetrended, so
    a drift that is present throughout the record is removed from the reference and
    left in the quantity being judged, and is refused.

    A stationary field exceeds the largest of N such windows with probability
    1/(N+1), which is what makes the retained record length, not a chosen number, set
    the rate at which the contract wrongly refuses.

    `biosphere/notes/equilibrium-trend-null.md` carries the measurement.
    """
    ncycle = policy["complete_forcing_cycles"]
    rate = float(policy["trend"]["cell_fraction"]["per_field_false_refusal_rate"])
    available = len(years) - window_years
    nwindow = available // window_years
    if nwindow < 1 or 1.0 / (nwindow + 1) > rate:
        needed = int(np.ceil(1.0 / rate) * window_years)
        raise EquilibriumWindowError(
            f"a retained record of {len(years)} years leaves {max(nwindow, 0)} "
            f"windows to measure the trending-cell null on, so a stationary field "
            f"would be refused at {1.0 / (max(nwindow, 0) + 1):.3f} against a "
            f"declared {rate:g}; retain at least {needed} years")
    start = available - nwindow * window_years
    block_years = years[start:available]
    if block_years != list(range(block_years[0], block_years[-1] + 1)):
        raise EquilibriumWindowError(
            "the record before the acceptance window has year gaps, so its "
            "windows cannot measure a null")
    block = cube[start:available]
    if not np.isfinite(block).all():
        raise EquilibriumWindowError(
            "the record before the acceptance window has missing cell-year rows")
    x = np.arange(block.shape[0], dtype=float)
    x -= x.mean()
    slope = np.einsum("t,tcf->cf", x, block) / float(np.sum(x * x))
    flat = block - x[:, None, None] * slope[None]
    windows = flat.reshape(nwindow, ncycle, cycle_years,
                           block.shape[1], block.shape[2]).mean(axis=2)
    fractions = np.stack([_cell_fraction(window, policy) for window in windows])
    limits = fractions.max(axis=0)
    return limits, {
        "statistic": policy["trend"]["cell_fraction"]["null_statistic"],
        "windows": int(nwindow),
        "first_year": int(block_years[0]), "last_year": int(block_years[-1]),
        "false_refusal_rate": 1.0 / (nwindow + 1),
        "per_field_false_refusal_rate": rate,
    }


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
    usable = (len(years) // cycle_years) * cycle_years
    if not np.isfinite(cube[:usable]).all():
        raise EquilibriumWindowError(
            f"{path}'s retained record has missing cell-year rows, so the memory "
            "time its acceptance depends on cannot be established")
    record = cube[:usable].reshape(
        usable // cycle_years, cycle_years, len(cells), len(names)).mean(axis=1)
    memory, unestablished = _memory_adequacy(record, names, policy)
    if unestablished:
        raise EquilibriumWindowError(
            f"{path}'s retained record is too short to judge: "
            f"{'; '.join(unestablished)}")
    limits, null = _cell_fraction_null(cube, years, window_years, cycle_years, policy)
    trends, rejected = _trend(cycle_means, policy, limits)
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
        "trend": {"rule": policy["trend"], "cell_fraction_null": null,
                  "memory": memory, "fields": trends, "verdict": "PASS"},
        "uncertainty": uncertainty,
        "peers": peer_records,
    }
    return ReducedTable(
        names=names,
        values={cell: mean[i] for i, cell in enumerate(cells)},
        temporal_std={cell: temporal_std[i] for i, cell in enumerate(cells)},
        report=report)
