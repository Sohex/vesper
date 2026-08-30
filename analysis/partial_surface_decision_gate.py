#!/usr/bin/env python3
"""Enforce SPAT-5's post-measurement representation decision."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "config/partial_surface.yaml"
REPORT = ROOT / "analysis/partial_surface_decision_gate_report.json"


def main() -> None:
    declaration = yaml.safe_load(CONTRACT.read_text(encoding="utf-8"))
    selection = declaration["selection"]
    measured_path = ROOT / selection["flux_measurement"]
    measured = json.loads(measured_path.read_text(encoding="utf-8"))
    decision = measured["selection"]
    convergence_path = ROOT / decision["reference"]
    convergence = json.loads(convergence_path.read_text(encoding="utf-8"))
    model_source = ROOT / "vendor/exoplasim/exoplasim/plasim/src"
    source = model_source / "oceanmod.f90"
    source_text = source.read_text(encoding="utf-8").lower()
    carrier_sources = {
        "builder": (ROOT / "exoplasim/scripts/build_boundary_conditions.py").read_text(),
        "runner": (ROOT / "exoplasim/scripts/run_exoplasim.py").read_text(),
        "restart": (ROOT / "exoplasim/scripts/restart_surface.py").read_text(),
        "surface": (ROOT / "vendor/exoplasim/exoplasim/plasim/src/surfmod.f90").read_text(),
        "state": (ROOT / "vendor/exoplasim/exoplasim/plasim/src/plasimmod.f90").read_text(),
    }
    checks = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})

    check("measurement contract",
          measured.get("contract_version") ==
              declaration["flux_bracket"]["contract_version"],
          str(measured.get("contract_version")))
    check("materiality test predates and matches the measurement",
          decision.get("materiality_test") == selection["materiality_test"],
          str(decision.get("materiality_test")))
    tolerance = float(
        convergence["criteria_provenance"]["storage_tolerance_w_m2"])
    check("materiality reference is unchanged",
          decision["state_storage_tolerance_w_m2"] == tolerance,
          f"{decision['state_storage_tolerance_w_m2']} vs {tolerance} W m-2")
    expected_tile = (
        decision["maximum_absolute_surface_energy_bracket_w_m2"] > tolerance)
    check("decision follows the declared inequality",
          decision["tile_model_selected"] == expected_tile
          and decision["model_form_bracket_selected"] == (not expected_tile),
          (f"bound {decision['maximum_absolute_surface_energy_bracket_w_m2']} "
           f"vs tolerance {tolerance} W m-2"))
    check("declaration records the measured decision",
          selection["tile_model_selected"] == decision["tile_model_selected"]
          and selection["model_form_bracket_selected"] ==
              decision["model_form_bracket_selected"],
          selection["tile_model_status"])

    tile = declaration.get("tile_model", {})
    check("selected tile operator is fully declared",
          tile.get("operator") == "one_atmosphere_column_two_surface_tiles"
          and tile.get("support", {}).get("land_fraction_source") ==
              "native_mesh_subaerial_area_fraction"
          and len(tile.get("state", {}).get("land", [])) == 4
          and len(tile.get("state", {}).get("ocean", [])) == 4
          and len(tile.get("exchange", {}).get("turbulent_fluxes", [])) == 4
          and len(tile.get("exchange", {}).get("radiative_fluxes", [])) == 3
          and len(tile.get("diagnostics", {}).get("retained_per_tile", [])) == 8
          and set(tile.get("acceptance", {})) == {
              "pure_land_reduction", "pure_ocean_reduction", "fraction_bounds",
              "extensive_closure", "storage_closure", "restart"},
          str(tile.get("operator")))

    boundary_report = json.loads((ROOT / "exoplasim/inputs/t21/"
                                  "boundary_conditions_report.json").read_text())
    fraction_path = (ROOT / "exoplasim/inputs/t21/"
                     "orogen_T21_surf_1720.sra")
    carrier_complete = (
        "LAND_FRACTION_CODE = 1720" in carrier_sources["builder"]
        and "LAND_FRACTION_SURFACE_CODES = {1720}" in carrier_sources["runner"]
        and '1720: "dlf"' in carrier_sources["restart"]
        and "call surfcode(1720,'dlf'" in carrier_sources["surface"]
        and "call surfcode(1720,'ylf'" in carrier_sources["surface"]
        and "call surfcode(1720,'xlf'" in carrier_sources["surface"]
        and "call mpsurfgp('dlf',dlf,NHOR,1)" in carrier_sources["surface"]
        and "call mpsurfgp('ylf',ylf,NHOR,1)" in
            (model_source / "oceanmod.f90").read_text()
        and "call mpsurfgp('xlf',xlf,NHOR,1)" in
            (model_source / "icemod.f90").read_text()
        and "real :: dlf(NHOR)" in carrier_sources["state"]
        and "real :: ylf(NHOR)" in
            (model_source / "oceanmod.f90").read_text()
        and "real :: xlf(NHOR)" in
            (model_source / "icemod.f90").read_text()
        and 1720 in boundary_report.get("codes", [])
        and boundary_report.get("land_fraction", {}).get("partial_cells", 0) > 0
        and fraction_path.is_file())
    check("selected fractional support reaches every model start",
          carrier_complete,
          (f"code 1720, {boundary_report.get('land_fraction', {}).get('partial_cells')} "
           "partial T21 cells; immutable boundary is re-read on restart"))

    hard_binary = ("where (yls(:) > 0.5)" in source_text
                   and "yls(:) = 1.0" in source_text
                   and "yls(:) = 0.0" in source_text)
    implementation_sources = {
        name: (model_source / name).read_text(encoding="utf-8").lower()
        for name in ("fluxmod.f90", "radmod.f90", "landmod.f90",
                     "seamod.f90", "outmod.f90")
    }
    required_seams = {
        "separate positive-fraction land state": (
            "spat5_tile_state" in implementation_sources["landmod.f90"]),
        "separate positive-fraction ocean state": (
            "spat5_tile_state" in implementation_sources["seamod.f90"]),
        "area-weighted turbulent exchange": (
            "spat5_tile_turbulent" in implementation_sources["fluxmod.f90"]),
        "area-weighted radiative exchange": (
            "spat5_tile_radiative" in implementation_sources["radmod.f90"]),
        "separate tile diagnostics": (
            "spat5_tile_diagnostics" in implementation_sources["outmod.f90"]),
        "tile restart records": all(
            "spat5_tile_restart" in implementation_sources[name]
            for name in ("landmod.f90", "seamod.f90")),
    }
    missing_seams = [name for name, present in required_seams.items()
                     if not present]
    check("selected representation is implemented",
          (not expected_tile or (not hard_binary and not missing_seams)),
          ("tile selected; "
           + ("oceanmod.f90 still hard-binarises the boundary; "
              if hard_binary else "")
           + ("missing " + ", ".join(missing_seams)
              if missing_seams else "all implementation seams present")))

    failed = [item for item in checks if not item["passed"]]
    report = {
        "contract_version": "vesper-partial-surface-decision-gate/1",
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "declaration": str(CONTRACT.relative_to(ROOT)),
        "measurement": str(measured_path.relative_to(ROOT)),
        "implementation": {
            "hard_binary_boundary": hard_binary,
            "required_seams": required_seams,
            "missing_seams": missing_seams,
        },
        "checks": checks,
        "summary": {"passed": len(checks) - len(failed), "failed": len(failed)},
        "verdict": "PASS" if not failed else "FAIL",
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    for item in checks:
        print(f"{'PASS' if item['passed'] else 'FAIL'}  {item['name']}: {item['detail']}")
    print(f"\n{report['verdict']}: {len(checks)-len(failed)} passed, "
          f"{len(failed)} failed; {REPORT.relative_to(ROOT)}")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
