#!/usr/bin/env python3
"""Validate the declared offline ocean-transport loop and write its gate report.

This is OCN-5's executable half. It does not run either model. It checks that
the decision in ``ocean/config/transport_loop.yaml`` still agrees with the
canonical pipeline graph and with the source seams on which the decision rests.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import yaml

from _paths import (ANALYSIS, CGENIE, LOOP_CONFIG, PIPELINE_CONFIG,
                    PLANET_CONFIG, PROJECT_ROOT)
from provenance import config_stamp

REPORT = ANALYSIS / "transport_loop_gate_report.json"
EXPECTED_STEPS = {
    "ocean_carbon_gate",
    "ocean_loop_gate",
    "ocean_flux_channel",
    "ocean_support",
    "ocean_forcing",
    "ocean_run",
    "ocean_convergence",
    "ocean_transport",
    "ocean_loop_convergence",
}
EXPECTED_CHAIN = [
    "baseline_climatology",
    "ocean_support",
    "ocean_forcing",
    "ocean_run",
    "ocean_convergence",
    "ocean_transport",
    "baseline_run",
    "baseline_climatology",
    "ocean_loop_convergence",
]


def compact(text: str) -> str:
    """Lowercase Fortran with comments and insignificant whitespace removed."""
    lines = [line.split("!", 1)[0] for line in text.splitlines()]
    return re.sub(r"\s+", "", "\n".join(lines)).lower()


def check(name: str, condition: bool, detail: str) -> dict:
    return {"check": name, "pass": bool(condition), "detail": detail}


def main() -> None:
    declaration = yaml.safe_load(LOOP_CONFIG.read_text(encoding="utf-8"))
    pipeline = yaml.safe_load(PIPELINE_CONFIG.read_text(encoding="utf-8"))
    planet = yaml.safe_load(PLANET_CONFIG.read_text(encoding="utf-8"))
    by_id = {step["id"]: step for step in pipeline["steps"]}
    loop_id = declaration.get("loop_id")

    genie = compact((CGENIE / "genie-main" / "genie.F").read_text(encoding="utf-8"))
    init = compact((CGENIE / "genie-main" / "initialise_genie.F").read_text(encoding="utf-8"))
    gold = compact((CGENIE / "genie-goldstein" / "src" / "fortran" /
                    "goldstein.F").read_text(encoding="utf-8"))
    gold_init = compact((CGENIE / "genie-goldstein" / "src" / "fortran" /
                         "initialise_goldstein.F").read_text(encoding="utf-8"))
    velc = compact((CGENIE / "genie-goldstein" / "src" / "fortran" /
                    "velc.f").read_text(encoding="utf-8"))

    architecture = declaration.get("architecture", {})
    procedure = declaration.get("forcing_procedure", {})
    criteria = declaration.get("exit", {}).get("criteria", [])
    criterion_ids = {row.get("id") for row in criteria}
    checks = [
        check("schema version", declaration.get("schema_version") == 1,
              "transport_loop schema version is 1"),
        check("loop is declared in the canonical graph", loop_id in pipeline.get("loops", {}),
              f"loop {loop_id!r} is present"),
        check("all ocean steps are registered", EXPECTED_STEPS <= set(by_id),
              f"missing: {sorted(EXPECTED_STEPS - set(by_id))}"),
        check("every ocean step belongs to the ocean loop",
              all(loop_id in ([by_id[s].get("loop")] if isinstance(by_id[s].get("loop"), str)
                              else by_id[s].get("loop", []))
                  for s in EXPECTED_STEPS if s in by_id),
              "each OCN-2/5/10/11 and execution row names loop O"),
        check("ocean physics is also inside loop A",
              all("A" in ([by_id[s].get("loop")] if isinstance(by_id[s].get("loop"), str)
                          else by_id[s].get("loop", []))
                  for s in EXPECTED_STEPS if s in by_id),
              "ocean transport can move the loop-A carve verdict"),
        check("offline full-flux architecture",
              architecture.get("coupling") == "offline_full_flux" and
              architecture.get("atmosphere_in_cgenie") == "none",
              "ExoPlaSim is the only atmosphere"),
        check("temperature restoring is absent",
              procedure.get("temperature_restoring") == "none",
              "the ocean is not restored toward the SST that forced it"),
        check("hosing is explicitly inert",
              procedure.get("freshwater_hosing") == {
                  "hosing_sv": 0.0,
                  "hosing_trend_sv_per_kyr": 0.0,
                  "years": 0,
              }, "all three always-reachable hosing controls are zero"),
        check("forcing gate exists in cGENIE",
              "if(flag_fluxatmos.and.(flag_ebatmos.or.flag_plasimatmos))then" in init and
              "if(flag_fluxatmos.and.flag_goldsteinocean)then" in genie,
              "the driven path is exclusive with both cGENIE atmospheres"),
        check("hosing remains on the source path",
              "callget_hosing(istep)" in gold and
              "namelist/ini_gold_nml/hosing,hosing_trend,nyears_hosing" in gold_init,
              "zero is a declared value, not an unreachable term"),
        check("solver rel is numerical under-relaxation",
              "rel*u1(1,i,j,k)+(1.0-rel)*u(1,i,j,k)" in velc and
              "rel*u1(2,i,j,k)+(1.0-rel)*u(2,i,j,k)" in velc,
              "rel mixes successive solver iterates and is not temperature restoring"),
        check("return set keeps ExoPlaSim sea ice authoritative",
              {row.get("name") for row in declaration.get("return_fields", [])} == {
                  "ocean_heat_flux_convergence",
                  "ocean_surface_eastward_velocity",
                  "ocean_surface_northward_velocity",
              }, "the ocean returns transport agents and no ice state"),
        check("every required exit axis is declared",
              {"heat_transport_change", "ocean_heat_integral",
               "sea_surface_temperature", "sea_ice_area",
               "atmospheric_energy_storage", "freshwater_closure",
               "salt_closure", "carve_drivers"} <= criterion_ids,
              f"criteria: {sorted(criterion_ids - {None})}"),
        check("transport loop has a finite refusal bound",
              isinstance(declaration.get("exit", {}).get("maximum_iterations"), int) and
              declaration["exit"]["maximum_iterations"] > 0 and
              declaration["exit"].get("on_maximum") ==
              "refuse_and_carry_an_unconverged_transport_bracket",
              "hitting the cap is not rounded into convergence"),
        check("provenance chain is chronological and complete",
              declaration.get("provenance_chain") == EXPECTED_CHAIN,
              " -> ".join(declaration.get("provenance_chain", []))),
        check("carbon feedback remains fail closed",
              architecture.get("carbon_feedback") ==
              "disabled_while_pco2_is_prescribed" and
              planet.get("atmosphere", {}).get("pCO2_bar") is not None,
              "marine carbon cannot alter a prescribed atmosphere"),
    ]

    report = {
        "source_build": planet.get("source_build"),
        "loop_id": loop_id,
        "checks": checks,
        "verdict": "PASS" if all(row["pass"] for row in checks) else "FAIL",
        **config_stamp(planet, "ocean/scripts/transport_loop_gate.py",
                       inputs=[LOOP_CONFIG, PIPELINE_CONFIG,
                               CGENIE / "genie-main" / "genie.F",
                               CGENIE / "genie-main" / "initialise_genie.F",
                               CGENIE / "genie-goldstein" / "src" / "fortran" /
                               "goldstein.F",
                               CGENIE / "genie-goldstein" / "src" / "fortran" /
                               "initialise_goldstein.F",
                               CGENIE / "genie-goldstein" / "src" / "fortran" /
                               "velc.f"]),
    }
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    for row in checks:
        print(f"[{' ok ' if row['pass'] else 'FAIL'}] {row['check']}: {row['detail']}")
    print(f"\n{report['verdict']}: {REPORT.relative_to(PROJECT_ROOT)}")
    raise SystemExit(0 if report["verdict"] == "PASS" else 1)


if __name__ == "__main__":
    main()
