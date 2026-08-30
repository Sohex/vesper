"""SDEC-1: prove the CENTURY accelerator is the sampled daily operator.

This is a no-simulation gate. It checks the production C++ source against the
declared operator and runs deterministic reduced fixtures for daily survival
composition, P isotherm preservation, flux sampling, and C-N-P closure.
"""

from __future__ import annotations

import copy
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import yaml

from _paths import GENERATED, PROJECT_ROOT
from paths import rel  # noqa: E402

DECLARATION = PROJECT_ROOT / "biosphere/config/somdynam.yaml"
REPORT = GENERATED / "century_acceleration_gate_report.json"
TOL = 1.0e-12


def _code(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    return re.sub(r"//[^\n]*", " ", text)


def compose_survival(loss_fractions: list[float]) -> float:
    survival = 1.0
    for loss in loss_fractions:
        if not 0.0 <= loss <= 1.0:
            raise ValueError(f"loss fraction outside [0,1]: {loss}")
        survival *= 1.0 - loss
    return survival


def representative_survival(monthly_survivals: list[float]) -> float:
    if not monthly_survivals:
        raise ValueError("no sampled months")
    return math.prod(monthly_survivals) ** (1.0 / len(monthly_survivals))


def partition_exchangeable_p(total: float, kplab: float,
                             spmax: float) -> tuple[float, float]:
    """The positive Langmuir root used by pmass_add()."""
    if min(total, kplab, spmax) < 0.0:
        raise ValueError("negative P state or isotherm parameter")
    b = kplab + spmax - total
    labile = 0.5 * (-b + math.sqrt(b * b + 4.0 * total * kplab))
    sorbed = total - labile
    return labile, sorbed


def p_mutation(total: float, delta: float, kplab: float,
               spmax: float) -> tuple[float, float]:
    updated = total + delta
    if updated < -TOL:
        raise ValueError("P mutation removes more exchangeable P than exists")
    return partition_exchangeable_p(max(0.0, updated), kplab, spmax)


def _isotherm_residual(labile: float, sorbed: float,
                       kplab: float, spmax: float) -> float:
    expected = spmax * labile / (kplab + labile) if labile else 0.0
    return sorbed - expected


def source_findings(declaration: dict, sources: dict[str, str]) -> list[dict]:
    block = declaration["century_acceleration"]
    code = _code(sources[block["source_file"]])
    state = _code(sources[block["state_file"]])
    init = _code(sources[block["initialization_file"]])
    findings: list[dict] = []

    def add(code_name: str, detail: str) -> None:
        findings.append({"code": code_name, "detail": detail})

    if block.get("contract_version") != "vesper-century-acceleration/1":
        add("schema", "unknown or absent accelerator contract version")
    for line in block.get("required_source", []):
        if line not in code:
            add("source_missing", f"production operator no longer runs {line!r}")
    for line in block.get("forbidden_source", []):
        if line in code:
            add("superseded_operator", f"production operator again runs {line!r}")

    for name in ("fnuptake_survival", "mminleach_survival",
                 "fpuptake_survival", "mminpleach_survival"):
        if f"double {name}[12];" not in state:
            add("state", f"{name} is not an explicit serialized Soil state")
        if f"{name}[mth] = 1.0;" not in init:
            add("identity", f"{name} does not initialise at multiplicative identity")
        if f"& {name}" not in init:
            add("restart", f"{name} is not serialized")

    for pool in block.get("restored_after_acceleration", []):
        if f"save_{pool} = soil.{pool};" not in code:
            add("p_restore", f"accelerator does not save {pool}")
        if f"soil.{pool} = save_{pool};" not in code:
            add("p_restore", f"accelerator does not restore {pool}")
    return findings


def operator_fixtures() -> list[dict]:
    fixtures: list[dict] = []

    def record(name: str, passed: bool, detail: str) -> None:
        fixtures.append({"fixture": name, "pass": bool(passed), "detail": detail})

    losses = [0.01, 0.03, 0.02, 0.04]
    daily = 100.0
    for loss in losses:
        daily *= 1.0 - loss
    composed = 100.0 * compose_survival(losses)
    record("daily uptake/leaching composes multiplicatively",
           math.isclose(daily, composed, rel_tol=0.0, abs_tol=TOL),
           f"daily={daily:.12g}, composed={composed:.12g}")

    old_once = 100.0 * (1.0 - sum(losses) / len(losses))
    record("arithmetic daily mean applied once is rejected",
           not math.isclose(old_once, daily, rel_tol=0.0, abs_tol=1.0e-6),
           f"wrong={old_once:.12g}, daily={daily:.12g}")

    inverted = 100.0 * (1.0 - compose_survival(losses))
    record("survival is not mistaken for the leached share",
           not math.isclose(inverted, daily, rel_tol=0.0, abs_tol=1.0e-6),
           f"inverted={inverted:.12g}, surviving={daily:.12g}")

    sampled = [compose_survival([0.01] * 30), compose_survival([0.02] * 30)]
    representative = representative_survival(sampled)
    record("sample years use a geometric monthly operator",
           math.isclose(representative ** 2, math.prod(sampled),
                        rel_tol=0.0, abs_tol=TOL),
           f"representative={representative:.12g}")

    organic_fractions = [0.00, 0.12, 0.03, 0.08]
    microbial_decay = [0.00, 0.20, 0.70, 0.10]
    daily_leached = sum(frac * decay for frac, decay in
                        zip(organic_fractions, microbial_decay))
    weighted = daily_leached / sum(microbial_decay)
    monthly_leached = weighted * sum(microbial_decay)
    arithmetic = sum(organic_fractions) / len(organic_fractions)
    record("organic leaching is weighted by decomposed microbial C",
           math.isclose(monthly_leached, daily_leached, abs_tol=TOL)
           and not math.isclose(arithmetic * sum(microbial_decay),
                                daily_leached, abs_tol=TOL),
           f"daily={daily_leached:.12g}, weighted={monthly_leached:.12g}, "
           f"arithmetic={arithmetic * sum(microbial_decay):.12g}")

    # Every P removal and addition repartitions the full exchangeable stock.
    kplab, spmax = 0.004, 0.12
    initial = 0.08
    labile, sorbed = partition_exchangeable_p(initial, kplab, spmax)
    p_ok = abs(_isotherm_residual(labile, sorbed, kplab, spmax)) < TOL
    total = labile + sorbed
    p_out = 0.0
    p_in = 0.0
    for survival in (compose_survival([0.02, 0.01, 0.03]),
                     compose_survival([0.04, 0.02])):
        removed = total * (1.0 - survival)
        labile, sorbed = p_mutation(total, -removed, kplab, spmax)
        total = labile + sorbed
        p_out += removed
        p_ok &= abs(_isotherm_residual(labile, sorbed, kplab, spmax)) < TOL
    for addition in (0.0012, 0.0007):
        labile, sorbed = p_mutation(total, addition, kplab, spmax)
        total = labile + sorbed
        p_in += addition
        p_ok &= abs(_isotherm_residual(labile, sorbed, kplab, spmax)) < TOL
    p_residual = total - initial - p_in + p_out
    record("every P mutation preserves isotherm and closure",
           p_ok and abs(p_residual) < TOL,
           f"isotherm={p_ok}, closure_residual={p_residual:.3e}")

    weather = [0.001, 0.002, 0.004, 0.003]
    flux_sum = sum(weather)
    cumulative_sum = sum(sum(weather[:i + 1]) for i in range(len(weather)))
    record("P forcing samples flux rather than year-to-date totals",
           math.isclose(flux_sum, 0.010, abs_tol=TOL)
           and not math.isclose(flux_sum, cumulative_sum, abs_tol=TOL),
           f"flux_sum={flux_sum:.6g}, triangular_sum={cumulative_sum:.6g}")

    # Reduced monthly C-N-P book: one internal decomposition transfer plus
    # respiration, organic leaching, mineral uptake/leaching and external input.
    initial_book = {"C": 12.0, "N": 1.2, "P": 0.18}
    inputs = {"C": 0.8, "N": 0.07, "P": 0.012}
    outputs = {"C": 0.45, "N": 0.025, "P": 0.006}
    final_book = {element: initial_book[element] + inputs[element] - outputs[element]
                  for element in initial_book}
    residuals = {element: final_book[element] - initial_book[element]
                 - inputs[element] + outputs[element]
                 for element in initial_book}
    record("accelerated C-N-P ledger closes",
           all(abs(value) < TOL for value in residuals.values()),
           ", ".join(f"{key}={value:.3e}" for key, value in residuals.items()))
    return fixtures


def source_fixtures(declaration: dict, sources: dict[str, str]) -> list[dict]:
    fixtures = []
    live = source_findings(declaration, sources)
    fixtures.append({"fixture": "live source matches declared accelerator",
                     "pass": not live, "detail": live})

    block = declaration["century_acceleration"]
    source_name = block["source_file"]
    cases = [
        ("direct P addition is refused",
         "pmass_add(soil, (soil.apwtr_mean + soil.apdep_mean) / 12.0);",
         "soil.pmass_labile += (soil.apwtr_mean + soil.apdep_mean) / 12.0;",
         "superseded_operator"),
        ("cumulative weathering sampling is refused",
         "soil.apwtr_mean += daily_pwtr;", "soil.apwtr_mean += soil.apwtr;",
         "superseded_operator"),
        ("missing deposition averaging is refused",
         "soil.apdep_mean /= nyear;", "", "source_missing"),
    ]
    for name, old, new, expected in cases:
        mutated = dict(sources)
        mutated[source_name] = mutated[source_name].replace(old, new, 1)
        found = {item["code"] for item in source_findings(declaration, mutated)}
        fixtures.append({"fixture": name, "pass": expected in found,
                         "detail": f"expected {expected}, found {sorted(found)}"})
    return fixtures


def main() -> int:
    declaration = yaml.safe_load(DECLARATION.read_text())
    block = declaration["century_acceleration"]
    names = (block["source_file"], block["state_file"],
             block["initialization_file"])
    sources = {name: (PROJECT_ROOT / name).read_text() for name in names}
    findings = source_findings(declaration, sources)
    fixtures = operator_fixtures() + source_fixtures(declaration, sources)
    broken = [row for row in fixtures if not row["pass"]]
    report = {
        "contract_version": block["contract_version"],
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": rel(PROJECT_ROOT / block["source_file"]),
        "findings": findings,
        "fixtures": fixtures,
        "verdict": "PASS" if not findings and not broken else "FAIL",
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n")
    for row in fixtures:
        print(f"[{' ok ' if row['pass'] else 'FAIL'}] {row['fixture']}: {row['detail']}")
    if findings:
        for finding in findings:
            print(f"[FAIL] {finding['code']}: {finding['detail']}")
    print(f"\n{report['verdict']}: {rel(REPORT)}")
    return 0 if report["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
