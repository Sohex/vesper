"""One equilibrium-window reducer for every LPJ-GUESS output consumer."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import struct
from typing import Iterable

import numpy as np
from scipy import stats
import yaml

from autocorrelation import (RELIABLE_SPAN_MULTIPLE, integrated_time,
                             mean_standard_error)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
POLICY_PATH = PROJECT_ROOT / "biosphere/config/equilibrium_window.yaml"
MEMORY_ESTIMATOR = "lib/autocorrelation.py:integrated_time"
DRIFT_ESTIMATOR = "lib/lpj_output.py:drift_bound"
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
    if policy.get("contract_version") != "vesper-lpj-equilibrium-window/5":
        raise EquilibriumWindowError("unsupported equilibrium-window contract")
    cycles = policy.get("complete_forcing_cycles")
    if not isinstance(cycles, int) or cycles < 3:
        raise EquilibriumWindowError("complete_forcing_cycles must be at least 3")
    trend = policy.get("trend", {})
    for key in ("relative_end_to_end_limit", "slope_standard_errors",
                "absolute_scale_floor"):
        if not isinstance(trend.get(key), (int, float)) or trend[key] <= 0:
            raise EquilibriumWindowError(f"trend.{key} must be positive")
    if policy.get("reported_span") != "whole_retained_record":
        raise EquilibriumWindowError(
            "reported_span names the span the reduced value is taken over, and "
            "this contract knows one: the whole retained record, which is the "
            "span it certifies")
    rate = trend.get("maximum_false_acceptance_rate")
    if not isinstance(rate, (int, float)) or not 0 < rate < 1:
        raise EquilibriumWindowError(
            "trend.maximum_false_acceptance_rate is the rate this contract "
            "controls and must lie in (0, 1)")
    if trend.get("drift_bound_estimator") != DRIFT_ESTIMATOR:
        raise EquilibriumWindowError(
            f"trend.drift_bound_estimator must name {DRIFT_ESTIMATOR}")
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
    _check_assessed(policy)
    return policy


def _check_assessed(policy: dict) -> None:
    """The assessed set is a list of consumer quantities and has to read as one.

    Every refusal here is a way the declaration could stop meaning what the
    contract claims: a quantity with no reader is not a consumer quantity, a
    quantity with no tolerance is the defect world-mxmr names, and a
    `relative_end_to_end_limit` that is not the tightest of them is a second
    number to keep in step with the first. The spin-up floor is derived at that
    tightest limit, so the two drifting apart would size a spin-up against a
    tolerance nothing is judged at.
    """
    assessed = policy.get("assessed")
    if not isinstance(assessed, dict):
        raise EquilibriumWindowError(
            "the contract needs an `assessed` block: the acceptance claim is "
            "about the quantities a consumer reads, and a column no consumer "
            "reads is a different claim from the one the tolerance was derived "
            "for")
    if assessed.get("reduction") != "sum":
        raise EquilibriumWindowError(
            "assessed.reduction names the one way a consumer forms a quantity "
            "out of columns, and this contract knows summation")
    quantities = assessed.get("quantities")
    if not isinstance(quantities, list) or not quantities:
        raise EquilibriumWindowError("assessed.quantities must list at least one")
    seen = set()
    for quantity in quantities:
        identity = quantity.get("id")
        if not isinstance(identity, str) or identity in seen:
            raise EquilibriumWindowError(
                f"every assessed quantity needs a unique id; {identity!r} is not")
        seen.add(identity)
        if not isinstance(quantity.get("table"), str):
            raise EquilibriumWindowError(f"{identity} names no table")
        explicit = quantity.get("columns")
        excluded = quantity.get("columns_excluding")
        if (explicit is None) == (excluded is None):
            raise EquilibriumWindowError(
                f"{identity} names its columns either explicitly or by exclusion, "
                "and exactly one of the two")
        for value in (explicit if explicit is not None else excluded):
            if not isinstance(value, str):
                raise EquilibriumWindowError(f"{identity} has a nonstring column")
        limit = quantity.get("relative_end_to_end_limit")
        if not isinstance(limit, (int, float)) or not 0 < limit < 1:
            raise EquilibriumWindowError(
                f"{identity} needs its own tolerance in (0, 1): applying one "
                "consumer's number to a quantity another consumer reads is the "
                "defect this block exists to end")
        readers = quantity.get("read_by")
        if not isinstance(readers, list) or not readers:
            raise EquilibriumWindowError(
                f"{identity} names no reader, so it is not a consumer quantity")
        if not isinstance(quantity.get("derivation"), str):
            raise EquilibriumWindowError(f"{identity} does not say where its "
                                         "tolerance came from")
    tightest = min(float(q["relative_end_to_end_limit"]) for q in quantities)
    declared = float(policy["trend"]["relative_end_to_end_limit"])
    if abs(tightest - declared) > 1.0e-12:
        raise EquilibriumWindowError(
            f"trend.relative_end_to_end_limit is {declared:g} and the tightest "
            f"assessed quantity's is {tightest:g}. It is not a separate number: "
            "it is what a spin-up has to leave below, and a spin-up precedes "
            "every quantity")


def assessed_quantities(policy: dict, table: str) -> list[dict]:
    """The consumer quantities this contract assesses on one output table."""
    return [quantity for quantity in policy["assessed"]["quantities"]
            if quantity["table"] == table]


def assess(cube: np.ndarray, names: list[str], quantities: list[dict]
           ) -> tuple[np.ndarray, list[str], np.ndarray]:
    """Replace a table's COLUMN axis with the consumer quantities built from it.

    `cube` is anything whose last axis is the table's columns; the same array
    comes back with that axis holding one entry per quantity, summed exactly as
    the consumer sums it. Every test in this contract then runs on the quantity
    a consumer will read rather than on a column that only enters one, which is
    what makes the acceptance claim and the consumer's claim the same claim.
    """
    cube = np.asarray(cube)
    if not quantities:
        return cube[..., :0], [], np.zeros(0)
    columns, labels, limits = [], [], []
    for quantity in quantities:
        explicit = quantity.get("columns")
        if explicit is not None:
            missing = [name for name in explicit if name not in names]
            if missing:
                raise EquilibriumWindowError(
                    f"{quantity['id']} needs columns {missing} that "
                    f"{quantity['table']} does not have")
            chosen = [names.index(name) for name in explicit]
        else:
            excluded = set(quantity["columns_excluding"])
            missing = [name for name in excluded if name not in names]
            if missing:
                raise EquilibriumWindowError(
                    f"{quantity['id']} excludes columns {missing} that "
                    f"{quantity['table']} does not have, so the exclusion no "
                    "longer describes the table it is taken on")
            chosen = [i for i, name in enumerate(names) if name not in excluded]
        if not chosen:
            raise EquilibriumWindowError(f"{quantity['id']} selects no column")
        columns.append(cube[..., chosen].sum(axis=-1))
        labels.append(quantity["id"])
        limits.append(float(quantity["relative_end_to_end_limit"]))
    return np.stack(columns, axis=-1), labels, np.asarray(limits, dtype=float)


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


def _cell_fraction(cycle_means: np.ndarray, policy: dict,
                   tolerances: np.ndarray | None = None) -> np.ndarray:
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
    if tolerances is None:
        tolerances = np.full(cycle_means.shape[2],
                             float(policy["trend"]["relative_end_to_end_limit"]))
    trending = ((relative > np.asarray(tolerances, dtype=float)[None, :])
                & (significance > float(policy["trend"]["slope_standard_errors"])))
    occupied = np.abs(intercept) > floor
    counts = occupied.sum(axis=0)
    return np.where(counts > 0, (trending & occupied).sum(axis=0)
                    / np.maximum(counts, 1), 0.0)


def _trend(cycle_means: np.ndarray, policy: dict,
           tolerances: np.ndarray,
           limits: np.ndarray) -> tuple[list[dict], bool]:
    """The PER-CELL half: drift that leaves the spatial mean flat.

    `_settled_within` owns the spatial mean and owns it over the whole record.
    What is left here is the one case a spatial mean cannot see, and it is why
    this half exists at all: opposed regional drifts that cancel. It is still a
    slope over the reported window against a limit measured from the run's own
    detrended windows, which means its size is calibrated rather than derived
    and its family false-refusal rate is well above the level it declares. That
    errs toward REFUSING, so a pass through it is evidence and a refusal through
    it is not; `biosphere/notes/equilibrium-trend-null.md` carries the numbers
    and the analytic replacement is tracked.
    """
    # cycle_means: cycle, cell, field
    nfield = cycle_means.shape[2]
    fractions = _cell_fraction(cycle_means, policy, tolerances)
    rejected = fractions > limits
    diagnostics = [{
        "trending_cell_fraction": float(fractions[i]),
        "trending_cell_fraction_limit": float(limits[i]),
        "cell_fraction_trending": bool(rejected[i]),
        "rejected": bool(rejected[i]),
    } for i in range(nfield)]
    return diagnostics, bool(rejected.any())


def drift_bound(series, alpha: float, floor: float) -> dict:
    """An upper confidence bound on a series' END-TO-END relative drift.

    WHAT THIS IS FOR, AND IT POINTS THE OPPOSITE WAY TO WHAT CAME BEFORE IT. An
    acceptance gate asserts that a run HAS SETTLED. A significance test that
    fails to reject "no drift" asserts nothing of the kind: it reports that it
    could not tell, and on a record shorter than its own memory time it can
    never tell about anything. The error such a gate must control is therefore
    letting a DRIFTING run through, and the instrument for that is an
    equivalence test -- bound the drift from above, pass only when the bound is
    inside the tolerance. Three properties follow from the direction alone, and
    each replaces a defect the earlier form could not repair.

    NO MULTIPLICITY CORRECTION IS NEEDED OR APPLIED. A run passes only when
    every assessed quantity passes, so by the intersection-union principle the
    run-level rate of accepting a field that truly drifts at the tolerance is
    bounded by `alpha` with nothing added. Contract 2 declared a per-field rate
    of 0.05 and measured a family rate of 0.34 over 64 fields, and expressing
    0.05 through its empirical null would have needed about 7100 retained
    cycles.

    NO SEPARATE MEMORY GUARD. A record too short to resolve the tolerance gives
    a wide bound, the bound exceeds the limit, and the field is refused by the
    test itself rather than by a span rule standing in front of it. Being a span
    rule is what made contract 3's guard self-referential.

    A CONVERGENT RECORD FLOOR, which `cycles_for_bound` inverts out of this.

    THE CONSTRUCTION. The record is split in half and the difference of the two
    half means is DOUBLED, because for a steady drift a half-to-half difference
    is half the end-to-end change; left unscaled it silently doubles the
    tolerance it is judged against. Each half mean carries the memory-corrected
    standard error `lib/autocorrelation.py` owns, at the memory time of the
    LARGER of the two halves' estimates -- a selection toward overestimates,
    which widens the bound and is the conservative direction here. Scatter and
    memory time are taken on the RAW half and not a detrended one, so a real
    drift inflates the error it is judged against rather than shrinking it.

    The bound is a t bound at Welch degrees of freedom built from each half's
    EFFECTIVE sample count, and that is where the span bar enters: a half
    carrying few effective samples earns a heavy critical value instead of a
    normal one. Putting a hard span bar on the whole record instead lets a half
    stand on five effective samples, which measured a family refusal rate of
    0.23 against a declared 0.05.
    """
    x = np.asarray(series, dtype=float)
    n = x.size
    m = n // 2
    if m < 3:
        raise EquilibriumWindowError(
            "a drift bound needs at least six cycles to halve")
    first, last = x[:m], x[n - m:]
    # THE UPPER END OF THE MEMORY TIME, NOT THE ESTIMATE. `integrated_time` is
    # biased low on a span carrying few independent samples -- against AR(1)
    # series of known memory time on 1253 samples it returns 86.1 for a true
    # 125.0 -- and a memory time read too small makes the standard error built
    # on it too small and the bound too narrow. On the point estimate the
    # declared acceptance rate held to a memory time of 175 cycles and reached
    # 0.109 at 400, which is why the upper end is used here. Taking the LARGER of
    # the two halves on top of that is the same conservatism applied to the split.
    # `biosphere/notes/equilibrium-trend-null.md` carries both sweeps.
    tau = max(integrated_time(first)["upper"], integrated_time(last)["upper"])
    variances = [float(np.var(half, ddof=1)) * tau / m for half in (first, last)]
    effective = m / tau
    per_half_df = max(effective - 1.0, 1.0)
    total = variances[0] + variances[1]
    if total > 0.0:
        df = total ** 2 / sum(v ** 2 / per_half_df for v in variances)
    else:
        df = per_half_df
    df = max(df, 1.0)
    drift = 2.0 * float(last.mean() - first.mean())
    standard_error = 2.0 * float(np.sqrt(total))
    scale = max(abs(float(x.mean())), float(floor))
    relative = abs(drift) / scale
    relative_error = standard_error / scale
    critical = float(stats.t.isf(alpha, df))
    return {"relative_drift": relative,
            "relative_standard_error": relative_error,
            "upper_bound": relative + critical * relative_error,
            "tau_cycles": float(tau),
            "effective_samples_per_half": float(effective),
            "degrees_of_freedom": float(df),
            "critical_value": critical,
            "half_cycles": int(m)}


def cycles_for_bound(relative_standard_error: float, cycles: int, tau: float,
                     alpha: float, limit: float, ceiling: float = 1.0e6) -> float:
    """The record a field needs before its drift bound can fall inside `limit`.

    THIS IS THE NUMBER THAT REPLACES THE SPAN MULTIPLE, and the reason it can is
    that it converges. `RELIABLE_SPAN_MULTIPLE * tau` asks for a record ten times
    a memory time read off the record itself, and the memory time this model
    reports grows with the window it is read on. Here the standard error falls as
    one over the root of the record while the memory time grows sublinearly with
    it, so the required length is reached rather than chased.

    The bound a SETTLED field of this scatter would show at length `n` is its
    expected absolute drift estimate plus the critical value times the standard
    error, and both terms scale with the same standard error:

        bound(n) = (E|z| + t(alpha, df(n))) * standard_error * sqrt(cycles / n)

    with `df(n) = n / tau - 2` for equal halves at a common memory time and
    `E|z| = sqrt(2 / pi)` for a standard normal. `bound` is decreasing in `n`, so
    the smallest `n` satisfying it is found by bisection on a doubled bracket.

    IT IS A FLOOR AND THE CALLER RECORDS IT AS ONE, for one reason: `tau` is held
    at its measured value while a longer record may read a larger one. That is a
    weaker self-reference than the span multiple's, because it moves the answer
    rather than preventing one.
    """
    if not np.isfinite(relative_standard_error) or relative_standard_error <= 0:
        return float(cycles)
    coefficient = float(relative_standard_error) * np.sqrt(float(cycles))
    expected_absolute = float(np.sqrt(2.0 / np.pi))

    def bound(n: float) -> float:
        df = max(n / float(tau) - 2.0, 1.0)
        critical = float(stats.t.isf(alpha, df))
        return (expected_absolute + critical) * coefficient / np.sqrt(n)

    if bound(float(cycles)) <= limit:
        return float(cycles)
    low, high = float(cycles), float(cycles) * 2.0
    while bound(high) > limit:
        high *= 2.0
        if high > ceiling:
            return float("inf")
    for _ in range(60):
        middle = 0.5 * (low + high)
        if bound(middle) <= limit:
            high = middle
        else:
            low = middle
    return high


def record_cycles_for_bound(series, cycles: int, alpha: float, limit: float,
                            floor: float) -> float:
    """The record a series of THIS scatter would need if it were settled.

    The question a refusal has to answer is how long the next run must be, and
    that is a question about a SETTLED field: a drift still present in the record
    inflates both the scatter and the memory time `drift_bound` reads off the raw
    halves, which is the right conservatism for the TEST and the wrong input for
    the LENGTH. So the standard error handed to `cycles_for_bound` is taken about
    the series' own linear fit, which removes the approach and leaves the
    variability the next record would still have to see through.
    """
    x = np.arange(np.asarray(series, dtype=float).size, dtype=float)
    flat = (np.asarray(series, dtype=float)
            - np.polyval(np.polyfit(x, series, 1), x) + np.mean(series))
    settled = drift_bound(flat, alpha, floor)
    return cycles_for_bound(settled["relative_standard_error"], cycles,
                            settled["tau_cycles"], alpha, limit)


def relaxation_time(series: np.ndarray, tau_memory: float,
                    blocks: int = 4) -> dict:
    """The e-folding time of an approach, measured WITHOUT its asymptote.

    WHY NOT A CURVE FIT. Fitting `a + b * exp(-t / tau)` needs the record to
    contain the turn-over: the asymptote `a` is a free parameter, and on a record
    shorter than the approach it lands outside the data and the fit says nothing.
    On this model's 1000-cycle record that happened for 34 of its 64 columns.

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
    exists to respect; the contraction must be RESOLVABLY below one, because a
    ratio whose own uncertainty reaches one is a record with no curvature in it
    and the timescale it implies is a lower bound rather than a value; and the
    second ratio must agree with the first inside a factor of two, because one
    exponential has one rate. A field failing any of them gets NO relaxation time
    and is not given a default.

    THE RESOLUTION CONDITION IS THE ONE THAT BITES HARDEST, and it is the same
    class of check as `drift_bound`'s: compare the instrument with the size of
    the effect before believing it. `tau = -Q / ln(ratio)` diverges as the ratio
    approaches one, so a ratio of 0.99 with an uncertainty of 0.35 spans several
    hundred cycles to unbounded and reads out as a confident five-figure
    timescale. On this model's simulated soil nitrogen it did exactly that, and a
    spin-up sized from it would have been thirty times what the same record
    supports.
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
    # The delta-method error on a ratio of two independent differences, which is
    # what decides whether the contraction is resolved at all.
    ratio_error = ratio * float(np.hypot(step_errors[0] / steps[0],
                                         step_errors[1] / steps[1]))
    if ratio + 2.0 * ratio_error >= 1.0:
        return {**base, "admissible": False, "ratio": ratio,
                "ratio_standard_error": ratio_error,
                "lower_bound_cycles": -q / float(np.log(
                    min(ratio + 2.0 * ratio_error, 1.0 - 1e-12)))
                if ratio + 2.0 * ratio_error < 1.0 else float("inf"),
                "reason": f"a contraction of {ratio:.3f} +/- {ratio_error:.3f} "
                          "reaches one, so this record carries no curvature and "
                          "the timescale it implies is a lower bound"}
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
    return {**base, "admissible": True, "ratio": ratio,
            "ratio_standard_error": ratio_error, "tau_cycles": tau,
            "second_estimate_cycles": second}


def timescale_report(run_dir: Path, *, policy_path: Path = POLICY_PATH) -> dict:
    """Both ecological timescales, per CONSUMER QUANTITY, for every table.

    Recorded on a run's acceptance artifact whatever its verdict, because a run
    that is refused for not having settled is exactly the run whose timescales say
    how long the next one has to be. `lib/run_lengths.py` reads this and states no
    number of its own, on the same terms as the climate relaxation bracket.

    `fields` HOLDS THE ASSESSED QUANTITIES AND NOTHING ELSE, because it is what
    sizes the next run: a record bought for a column no consumer reads is a
    record bought for a claim nobody makes. Every column's own timescales are in
    `columns` beside it, which nothing reads and a reader can.
    """
    run_dir = Path(run_dir)
    policy = read_policy(policy_path)
    floor = float(policy["trend"]["absolute_scale_floor"])
    tightest = float(policy["trend"]["relative_end_to_end_limit"])
    alpha = float(policy["trend"]["maximum_false_acceptance_rate"])
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
        quantities = assessed_quantities(policy, output)
        assessed, labels, limits = (
            assess(record, names, quantities) if quantities
            else (record[:, :, :0], [], np.zeros(0)))

        def timescales(spatial: np.ndarray, limit: float) -> dict | None:
            if (abs(float(spatial.mean())) <= floor
                    or not np.isfinite(spatial).all() or spatial.std() == 0):
                return None
            flat = (spatial - np.polyval(np.polyfit(x, spatial, 1), x)
                    + spatial.mean())
            memory = integrated_time(flat)
            bound = drift_bound(spatial, alpha, floor)
            return {
                "memory": {"tau_cycles": float(memory["tau"]),
                           "effective_samples": float(memory["effective_sample_size"]),
                           "reliable": bool(memory["reliable"]),
                           "lag1": float(memory["lag1"])},
                "relaxation": relaxation_time(spatial, memory["tau"]),
                # The RECORD this quantity needs, recorded on a refusal as well
                # as on a pass because a run refused for not resolving its own
                # drift is exactly the run that says how long the next one has
                # to be. `lib/run_lengths.py` reads it and states no number of
                # its own.
                "drift": {**bound,
                          "relative_end_to_end_limit": float(limit),
                          "settled": bool(bound["upper_bound"] <= limit),
                          "record_cycles_for_bound": float(
                              record_cycles_for_bound(spatial, span, alpha,
                                                      limit, floor))},
            }

        fields = {}
        for index, (label, limit) in enumerate(zip(labels, limits)):
            entry = timescales(assessed[:, :, index].mean(axis=1), float(limit))
            if entry is not None:
                entry["read_by"] = list(quantities[index]["read_by"])
                fields[label] = entry
        columns = {}
        for index, name in enumerate(names):
            entry = timescales(record[:, :, index].mean(axis=1), tightest)
            if entry is not None:
                columns[name] = entry
        tables[output] = {"record_cycles": int(span),
                          "forcing_cycle_years": int(cycle_years),
                          "fields": fields, "columns": columns}
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


def _bound_series(spatial: np.ndarray, span: int, limit: float, alpha: float,
                  floor: float) -> tuple[dict | None, float]:
    """The drift bound on one series, and the record a settled one would need."""
    level = float(abs(spatial.mean()))
    if level <= floor or not np.isfinite(spatial).all() or spatial.std() == 0:
        return None, float("nan")
    bound = drift_bound(spatial, alpha, floor)
    needed = record_cycles_for_bound(spatial, span, alpha, limit, floor)
    # What the REPORTED value is worth, in the units the tolerance is in, so a
    # consumer never has to reconstruct it. This is the one standard error of a
    # mean over a series with memory, at the same memory time the bound was taken
    # at, over the same span the bound certifies.
    bound["reported_mean_relative_standard_error"] = float(
        mean_standard_error(spatial, bound["tau_cycles"]) / level)
    bound["record_cycles"] = int(span)
    bound["record_cycles_for_bound"] = float(needed)
    return bound, needed


def _settled_within(record: np.ndarray, names: list[str], policy: dict,
                    quantities: list[dict]) -> tuple[list[dict], list[dict], list[str]]:
    """Is every CONSUMER QUANTITY's drift demonstrably inside its own tolerance?

    THE GLOBAL HALF OF THE CONTRACT, and it runs on the WHOLE RETAINED RECORD
    rather than on the reported window. Those are two different spans and only
    one number used to name both: the window is how much of the run a consumer
    reads a value over, and the record is how much of it the equilibrium claim
    is made from. A drift test on the window was a slope fitted inside one
    memory time, which is what `drift_bound` replaces.

    IT RUNS ON THE QUANTITY AND NOT ON THE COLUMN. A consumer sums columns and
    reads the sum, and the sum's own series carries each column's drift at that
    column's share, with every cancellation and every reinforcement in it. So the
    bound applied to the sum bounds the sum exactly and the columns beneath it
    owe nothing further; bounding each separately at the sum's tolerance is a
    stronger claim than any consumer makes, and it cost this contract a
    twenty-seven-fold record. The columns keep a bound as a DIAGNOSTIC, returned
    separately, so what is outside the claim is visible rather than absent.

    A quantity is settled when the upper bound on its end-to-end relative drift
    is inside ITS OWN tolerance, and it is REFUSED otherwise, whether the bound
    is wide because the quantity is drifting or wide because the record is short.
    Those two are the same refusal on purpose: an acceptance gate that passes a
    run it cannot resolve is asserting something it did not measure. The refusal
    names which of the two it is, by carrying both the drift estimate and the
    record that would answer it.
    """
    floor = float(policy["trend"]["absolute_scale_floor"])
    alpha = float(policy["trend"]["maximum_false_acceptance_rate"])
    span = record.shape[0]
    assessed, labels, limits = assess(record, names, quantities)
    diagnostics, refused = [], []
    for index, (label, limit) in enumerate(zip(labels, limits)):
        spatial = assessed[:, :, index].mean(axis=1)
        bound, needed = _bound_series(spatial, span, float(limit), alpha, floor)
        if bound is None:
            diagnostics.append({
                "field": label, "assessed": False,
                "relative_end_to_end_limit": float(limit),
                "reason": "the spatial mean is zero or exactly constant, so it "
                          "carries no drift to bound"})
            continue
        settled = bool(bound["upper_bound"] <= float(limit))
        diagnostics.append({"field": label, "assessed": True, "settled": settled,
                            "relative_end_to_end_limit": float(limit),
                            "read_by": list(quantities[index]["read_by"]),
                            "derivation": quantities[index]["derivation"],
                            **bound})
        if settled:
            continue
        length = ("no finite record" if not np.isfinite(needed)
                  else f"{needed:.0f} cycles")
        refused.append(
            f"{label} (drift {bound['relative_drift']:.4f} +/- "
            f"{bound['relative_standard_error']:.4f} bounds at "
            f"{bound['upper_bound']:.4f} against a limit of {float(limit):g}; "
            f"memory time {bound['tau_cycles']:.1f} cycles leaves "
            f"{bound['effective_samples_per_half']:.1f} effective samples per "
            f"half of a {span}-cycle record, and a SETTLED quantity of this "
            f"scatter would resolve the limit at {length})")
    # THE COLUMNS, at the tightest limit any assessed quantity carries, so the
    # diagnostic is comparable across tables. Nothing is gated on these: they
    # exist so a reader can see what the claim does NOT cover.
    tightest = float(policy["trend"]["relative_end_to_end_limit"])
    columns = []
    for index, name in enumerate(names):
        spatial = record[:, :, index].mean(axis=1)
        bound, _ = _bound_series(spatial, span, tightest, alpha, floor)
        if bound is None:
            columns.append({"field": name, "bounded": False,
                            "reason": "the spatial mean is zero or exactly "
                                      "constant, so it carries no drift to bound"})
            continue
        columns.append({"field": name, "bounded": True,
                        "compared_against": tightest, **bound})
    return diagnostics, columns, refused


def _cell_fraction_null(cube: np.ndarray, years: list[int], window_years: int,
                        cycle_years: int, policy: dict,
                        tolerances: np.ndarray | None = None
                        ) -> tuple[np.ndarray, dict]:
    """Per quantity, the trending-cell fraction this run reaches with no trend left.

    The limit on the trending-cell fraction cannot be one declared number. Measured
    on a run at fixed forcing, the fraction a field reaches when nothing is drifting
    spans three orders of magnitude across a table's columns, because it is set by
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
    fractions = np.stack([_cell_fraction(window, policy, tolerances)
                          for window in windows])
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
    # THE REPORTED VALUE IS TAKEN OVER THE SPAN THE CONTRACT CERTIFIES, and that
    # is the whole retained record. Reporting a ten-cycle mean while certifying
    # the record hands a consumer a number whose own sampling error is larger
    # than the drift the certificate refuses: the memory time of most assessed
    # fields is tens to hundreds of cycles, so a ten-cycle mean is one effective
    # sample and cannot be known better than the field's marginal scatter, which
    # for ten of this model's columns exceeds `relative_end_to_end_limit` outright.
    # The window keeps one job, and it is the per-cell half's, whose empirical
    # null is built from windows and needs many of them.
    reported = cube[:usable]
    mean = reported.mean(axis=0)
    temporal_std = reported.std(axis=0, ddof=1)
    # THE CLAIM IS ABOUT THE QUANTITIES A CONSUMER READS, so both halves run on
    # the sums the consumers form rather than on this table's columns. A table
    # nothing reads carries no quantity and makes no claim, which is a pass and
    # is recorded as one rather than as an assessment of nothing.
    quantities = assessed_quantities(policy, path.name)
    drift, column_drift, unsettled = _settled_within(
        record, names, policy, quantities)
    if unsettled:
        raise EquilibriumWindowError(
            f"{path}'s retained record does not bound the drift of the "
            f"quantities its consumers read inside their own tolerances: "
            f"{'; '.join(unsettled)}")
    trends, null = [], None
    if quantities:
        cycle_quantities, labels, tolerances = assess(
            cycle_means, names, quantities)
        null_cube, _, _ = assess(cube, names, quantities)
        limits, null = _cell_fraction_null(
            null_cube, years, window_years, cycle_years, policy, tolerances)
        trends, rejected = _trend(cycle_quantities, policy, tolerances, limits)
        for label, diagnostic in zip(labels, trends):
            diagnostic["field"] = label
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
        "policy": str(policy_path.relative_to(PROJECT_ROOT)
                      if policy_path.is_relative_to(PROJECT_ROOT)
                      else policy_path),
        "policy_sha256": sha256(policy_path),
        "identity": {"table": str(path), "run_manifest": str(manifest_path),
                     "run_id": manifest.get("run_id"),
                     "root_seed": manifest.get("stochastic_randomness", {}).get(
                         "root_seed", manifest.get("physical", {}).get("root_seed")),
                     "npatch": manifest.get("physical", {}).get("npatch")},
        "forcing": forcing_source,
        # The span `values` and `temporal_std` are taken over, which is the span
        # the drift bound certifies. A consumer reads this one.
        "reported": {"first_year": years[0], "last_year": years[usable - 1],
                     "annual_values": int(usable),
                     "complete_forcing_cycles": int(usable // cycle_years),
                     "forcing_cycle_years": cycle_years},
        # The per-cell half's window, and nothing else's.
        "window": {"first_year": selected_years[0], "last_year": selected_years[-1],
                   "annual_values": window_years,
                   "complete_forcing_cycles": policy["complete_forcing_cycles"],
                   "forcing_cycle_years": cycle_years},
        # WHAT THE PASS ASSERTS, in the words the claim is made in: these
        # quantities, each inside the tolerance its own consumer owes. It is not
        # an assertion that the simulated biosphere has settled, and
        # `column_drift` carries every column's bound so a reader can see what
        # the claim leaves out rather than having to notice its absence.
        "trend": {"rule": policy["trend"], "cell_fraction_null": null,
                  "record_cycles": int(record.shape[0]), "drift": drift,
                  "assessed": [q["id"] for q in quantities],
                  "claim": ("the quantities named in `drift` are inside the "
                            "tolerance each owes its own consumer over the whole "
                            "retained record; no claim is made about a column "
                            "no consumer reads"),
                  "column_drift": column_drift,
                  "fields": trends, "verdict": "PASS"},
        "uncertainty": uncertainty,
        "peers": peer_records,
    }
    return ReducedTable(
        names=names,
        values={cell: mean[i] for i, cell in enumerate(cells)},
        temporal_std={cell: temporal_std[i] for i, cell in enumerate(cells)},
        report=report)
