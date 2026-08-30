#!/usr/bin/env python3
"""Keep marine carbon feedback closed while atmospheric CO2 is prescribed.

OCN-16 does not activate a carbon cycle. This gate proves the opposite: every
ExoPlaSim entry point explicitly disables its resident CO2 evolution, the ocean
transport return carries no carbon term, and changing the prescribed premise
without first supplying a conservative replacement contract is a refusal.
"""

from __future__ import annotations

import ast
import copy
import json
import math

import yaml

from _paths import (ANALYSIS, CARBON_CONFIG, LOOP_CONFIG, PIPELINE_CONFIG,
                    PLANET_CONFIG, PROJECT_ROOT)
from provenance import config_stamp

REPORT = ANALYSIS / "carbon_feedback_gate_report.json"
RUN_SOURCES = [
    PROJECT_ROOT / "exoplasim" / "scripts" / "run_exoplasim.py",
    PROJECT_ROOT / "exoplasim" / "scripts" / "continue_exoplasim.py",
    PROJECT_ROOT / "exoplasim" / "scripts" / "run_stellar_cycle.py",
]
REQUIRED_REOPEN = {"anut-6", "ocn-13", "ocn-14", "volc-8"}
CARBON_RETURN_WORDS = {
    "co2", "carbon", "dic", "alkalinity", "carbonate", "gas_exchange",
}


def check(name: str, condition: bool, detail: str) -> dict:
    return {"check": name, "pass": bool(condition), "detail": detail}


def configure_calls_are_fixed(source: str) -> tuple[bool, str]:
    """Return whether every model.configure call fixes both carbon switches."""
    tree = ast.parse(source)
    calls = [node for node in ast.walk(tree)
             if isinstance(node, ast.Call)
             and isinstance(node.func, ast.Attribute)
             and node.func.attr == "configure"]
    if not calls:
        return False, "no configure call found"
    for call in calls:
        keywords = {row.arg: row.value for row in call.keywords if row.arg}
        if "pCO2" not in keywords:
            return False, "configure call omits prescribed pCO2"
        for key in ("co2weathering", "evolveco2"):
            value = keywords.get(key)
            if not (isinstance(value, ast.Constant) and value.value is False):
                return False, f"configure call does not set {key}=False"
    return True, f"{len(calls)} configure call(s) fix both switches false"


def evaluate(carbon: dict, planet: dict, transport: dict, pipeline: dict,
             source_guards: dict[str, tuple[bool, str]]) -> list[dict]:
    pco2 = planet.get("atmosphere", {}).get("pCO2_bar")
    active = carbon.get("active_boundary", {})
    reopening = carbon.get("reopening_contract", {})
    conservation = reopening.get("conservation", {})
    returns = [str(row.get("name", "")).lower()
               for row in transport.get("return_fields", [])]
    carbon_returns = sorted(name for name in returns
                            if any(word in name for word in CARBON_RETURN_WORDS))
    steps = {row.get("id"): row for row in pipeline.get("steps", [])}
    loop_needs = steps.get("ocean_loop_gate", {}).get("needs", [])
    return [
        check("schema and active mode",
              carbon.get("schema_version") == 1 and
              carbon.get("mode") == "fail_closed_prescribed_atmospheric_pco2",
              str(carbon.get("mode"))),
        check("prescribed atmospheric pCO2 exists",
              isinstance(pco2, (int, float)) and math.isfinite(float(pco2))
              and float(pco2) > 0.0,
              f"pCO2_bar={pco2!r}"),
        check("ocean cannot modify the prescribed atmosphere",
              active.get("ocean_air_sea_co2_exchange") == "disabled" and
              active.get("ocean_may_modify_atmospheric_co2") is False and
              active.get("atmospheric_carbon_reservoir") ==
              "absent_while_prescribed",
              "no active atmospheric carbon reservoir or exchange"),
        check("transport architecture agrees",
              transport.get("architecture", {}).get("carbon_feedback") ==
              "disabled_while_pco2_is_prescribed",
              str(transport.get("architecture", {}).get("carbon_feedback"))),
        check("ocean return contains no carbon term", not carbon_returns,
              f"forbidden returns: {carbon_returns}"),
        check("all ExoPlaSim entry points pin carbon evolution off",
              all(value[0] for value in source_guards.values()),
              "; ".join(f"{name}: {value[1]}"
                         for name, value in source_guards.items())),
        check("carbon gate is on the ocean-loop path",
              "ocean_carbon_gate" in steps and
              "ocean_carbon_gate" in loop_needs,
              f"ocean_loop_gate needs={loop_needs}"),
        check("reopening dependencies are complete",
              set(reopening.get("required_closed_issues", [])) == REQUIRED_REOPEN,
              str(reopening.get("required_closed_issues", []))),
        check("reopening requires permission, a loop, and convergence",
              reopening.get("requires_explicit_permission") is True and
              reopening.get("requires_new_pipeline_loop") is True and
              reopening.get("requires_declared_convergence_predicate") is True,
              "all three escalation gates must be true"),
        check("air-sea transfer is conservative",
              conservation.get("transfer_identity") ==
              "atmospheric_debit_mol_c + ocean_credit_mol_c == 0" and
              conservation.get("prohibit_unbalanced_ocean_uptake") is True and
              conservation.get("require_restart_continuity") is True,
              str(conservation.get("transfer_identity"))),
    ]


def mutation_fixtures(carbon: dict, planet: dict, transport: dict, pipeline: dict,
                      guards: dict[str, tuple[bool, str]]) -> list[dict]:
    fixtures = []

    def expect_refusal(name: str, c=carbon, p=planet, t=transport, g=guards) -> None:
        passed = all(row["pass"] for row in evaluate(c, p, t, pipeline, g))
        fixtures.append({"fixture": name, "pass": not passed})

    changed = copy.deepcopy(planet)
    changed["atmosphere"]["pCO2_bar"] = None
    expect_refusal("missing prescribed pCO2", p=changed)
    changed = copy.deepcopy(carbon)
    changed["mode"] = "active_air_sea_exchange"
    expect_refusal("undeclared active mode", c=changed)
    changed = copy.deepcopy(transport)
    changed["architecture"]["carbon_feedback"] = "enabled"
    expect_refusal("transport enables carbon", t=changed)
    changed = copy.deepcopy(transport)
    changed["return_fields"].append({"name": "atmospheric_co2_tendency"})
    expect_refusal("unbalanced carbon return", t=changed)
    changed = copy.deepcopy(guards)
    changed[next(iter(changed))] = (False, "fixture removed evolveco2=False")
    expect_refusal("run entry point evolves CO2", g=changed)
    changed = copy.deepcopy(carbon)
    changed["reopening_contract"]["conservation"][
        "prohibit_unbalanced_ocean_uptake"] = False
    expect_refusal("unbalanced uptake permitted", c=changed)
    changed = copy.deepcopy(carbon)
    changed["reopening_contract"]["required_closed_issues"].remove("anut-6")
    expect_refusal("river export dependency omitted", c=changed)
    return fixtures


def main() -> None:
    carbon = yaml.safe_load(CARBON_CONFIG.read_text(encoding="utf-8"))
    planet = yaml.safe_load(PLANET_CONFIG.read_text(encoding="utf-8"))
    transport = yaml.safe_load(LOOP_CONFIG.read_text(encoding="utf-8"))
    pipeline = yaml.safe_load(PIPELINE_CONFIG.read_text(encoding="utf-8"))
    guards = {path.relative_to(PROJECT_ROOT).as_posix():
              configure_calls_are_fixed(path.read_text(encoding="utf-8"))
              for path in RUN_SOURCES}
    checks = evaluate(carbon, planet, transport, pipeline, guards)
    fixtures = mutation_fixtures(carbon, planet, transport, pipeline, guards)
    passed = all(row["pass"] for row in checks + fixtures)
    report = {
        "source_build": planet.get("source_build"),
        "mode": carbon.get("mode"),
        "checks": checks,
        "mutation_fixtures": fixtures,
        "verdict": "PASS" if passed else "FAIL",
        **config_stamp(planet, "ocean/scripts/carbon_feedback_gate.py",
                       inputs=[CARBON_CONFIG, PLANET_CONFIG, LOOP_CONFIG,
                               PIPELINE_CONFIG, *RUN_SOURCES]),
    }
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    for row in checks:
        print(f"[{' ok ' if row['pass'] else 'FAIL'}] {row['check']}: {row['detail']}")
    for row in fixtures:
        print(f"[{' ok ' if row['pass'] else 'FAIL'}] fixture: {row['fixture']}")
    print(f"\n{report['verdict']}: {REPORT.relative_to(PROJECT_ROOT)}")
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
