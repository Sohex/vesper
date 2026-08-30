#!/usr/bin/env python3
"""Exercise SPAT-1's spatial-support contract against positive/negative cases."""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))

from spatial_support import SpatialContractError, declaration, validate_contract  # noqa: E402

REPORT = ROOT / "analysis" / "spatial_support_gate_report.json"
H0 = "0" * 64
H1 = "1" * 64
EXPECTED_SUPPORT_KINDS = {
    "terrain_native_mesh",
    "atmosphere_gaussian_grid",
    "ocean_cgenie_grid",
    "hydrologic_unit",
    "soil_response_unit",
    "land_ecological_response_unit",
    "marine_ecological_response_unit",
    "lpj_grid_cell",
    "lpj_demographic_patch",
    "comparison_support",
}


def reference_contract() -> dict:
    support = {
        "id": "atmosphere_t42",
        "kind": "atmosphere_gaussian_grid",
        "dimensions": {"latitude": 64, "longitude": 128},
        "geometry_sha256": H0,
        "coordinates_sha256": H1,
        "native_measures": [
            {"kind": "area", "variable": "cell_area", "units": "m2"},
        ],
    }
    aggregation = {
        "operator": "identity",
        "operator_version": "lib.gridding/1",
        "source_support": "atmosphere_t42",
        "destination_support": "atmosphere_t42",
        "nonlinear_order": "not_applicable",
    }
    return {
        "contract_version": "vesper-spatial-support/1",
        "artifact_id": "synthetic/atmosphere-grid",
        "supports": [support],
        "fields": [
            {
                "name": "land_fraction",
                "semantics": "dimensionless_fraction",
                "units": "1",
                "support_id": "atmosphere_t42",
                "spatial_dimensions": ["latitude", "longitude"],
                "effective_measure": {"mode": "native", "measure_kind": "area"},
                "aggregation": aggregation,
            },
            {
                "name": "soil_carbon",
                "semantics": "extensive_total",
                "units": "kg C",
                "support_id": "atmosphere_t42",
                "spatial_dimensions": ["latitude", "longitude"],
                "effective_measure": {
                    "mode": "native_times_fraction",
                    "measure_kind": "area",
                    "fraction_field": "land_fraction",
                },
                "aggregation": aggregation,
            },
        ],
        "time_support": {
            "semantics": "interval_mean",
            "calendar": "vesper_orbit",
            "absolute_time_units": "seconds since model epoch",
            "bounds_variable": "time_bounds",
            "interval_closure": "[start,end)",
        },
        "provenance": {
            "source_artifacts": ["synthetic://mesh"],
            "operator_versions": {"identity": "lib.gridding/1"},
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=REPORT)
    args = parser.parse_args()
    vocab = declaration()
    base = reference_contract()
    checks: list[dict] = [
        {
            "check": "schema and contract version are fixed",
            "pass": vocab.get("schema_version") == 1 and
                    vocab.get("contract_version") == "vesper-spatial-support/1",
        },
        {
            "check": "every distinct project support is named",
            "pass": set(vocab.get("support_kinds", {})) == EXPECTED_SUPPORT_KINDS,
            "declared": sorted(vocab.get("support_kinds", {})),
        },
        {
            "check": "LPJ patches cannot masquerade as response units",
            "pass": vocab["support_kinds"]["lpj_demographic_patch"].get(
                "spatially_representative") is False,
        },
    ]

    def positive(name: str, contract: dict) -> None:
        try:
            identity = validate_contract(contract, vocab)
            checks.append({"check": name, "pass": True, "identity": identity})
        except SpatialContractError as exc:
            checks.append({"check": name, "pass": False, "detail": str(exc)})

    def negative(name: str, mutate) -> None:
        candidate = copy.deepcopy(base)
        mutate(candidate)
        try:
            validate_contract(candidate, vocab)
        except SpatialContractError as exc:
            checks.append({"check": name, "pass": True, "refusal": str(exc)})
        else:
            checks.append({"check": name, "pass": False,
                           "detail": "malformed contract was accepted"})

    positive("complete chronological contract", base)
    negative("shape is not support identity",
             lambda c: c["supports"][0].pop("coordinates_sha256"))
    negative("unknown support kind refuses",
             lambda c: c["supports"][0].update(kind="generic_resolution"))
    negative("non-positive dimensions refuse",
             lambda c: c["supports"][0]["dimensions"].update(latitude=0))
    negative("native measure is mandatory",
             lambda c: c["supports"][0].pop("native_measures"))
    negative("field shape names its support dimensions",
             lambda c: c["fields"][0].update(spatial_dimensions=["cell"]))
    negative("field semantics are closed",
             lambda c: c["fields"][1].update(semantics="number"))
    negative("effective fraction must exist",
             lambda c: c["fields"][1]["effective_measure"].update(
                 fraction_field="rootable_fraction"))
    negative("fraction must share support",
             lambda c: c["fields"][0].update(support_id="missing"))
    negative("operator version is mandatory",
             lambda c: c["fields"][0]["aggregation"].pop("operator_version"))
    negative("chronology needs absolute time",
             lambda c: c["time_support"].pop("absolute_time_units"))
    negative("interval closure is explicit",
             lambda c: c["time_support"].update(interval_closure="ambiguous"))
    negative("demographic patches are not spatial units", lambda c: (
        c["supports"][0].update(kind="lpj_demographic_patch"),
        c["supports"][0].pop("spatially_representative", None)))
    negative("vector crossing declares rotation", lambda c: (
        c["fields"][0].update(semantics="vector_component"),
        c["fields"][0]["aggregation"].update(operator="covered_area_mean")))
    negative("vector basis is explicit", lambda c: (
        c["fields"][0].update(semantics="vector_component"),
        c["fields"][0]["aggregation"].update(operator="identity")))

    changed = copy.deepcopy(base)
    changed["supports"][0]["coordinates_sha256"] = "2" * 64
    identity_changes = validate_contract(base, vocab) != validate_contract(changed, vocab)
    checks.append({"check": "coordinate identity changes the contract digest",
                   "pass": identity_changes})

    report = {
        "schema_version": vocab.get("schema_version"),
        "contract_version": vocab.get("contract_version"),
        "checks": checks,
        "verdict": "PASS" if all(row["pass"] for row in checks) else "FAIL",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    for row in checks:
        print(f"[{' ok ' if row['pass'] else 'FAIL'}] {row['check']}")
    print(f"\n{report['verdict']}: {args.output.relative_to(ROOT)}")
    raise SystemExit(0 if report["verdict"] == "PASS" else 1)


if __name__ == "__main__":
    main()
