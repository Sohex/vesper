#!/usr/bin/env python3
"""Run SPAT-10's assessed-conversion schema and no-simulation fixtures."""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))

import gridding  # noqa: E402
import remap  # noqa: E402
from spatial_support import (SpatialContractError, canonical_digest,  # noqa: E402
                             conversion_declaration,
                             grid_support_contract,
                             validate_contract,
                             validate_conversion_assessment)

REPORT = ROOT / "analysis" / "spatial_conversion_gate_report.json"


def support_contract(spec: gridding.GridSpec, sid: str, kind: str) -> dict:
    measures = [{"kind": "area", "variable": "cell_area", "units": "m2"}]
    if kind == "ocean_cgenie_grid":
        measures.append({"kind": "volume", "variable": "cell_volume", "units": "m3"})
    return grid_support_contract(
        spec, f"synthetic/{sid}", sid, kind,
        source_artifacts=["synthetic://fixture"], native_measures=measures)


def pass_row(residual: float, tolerance: float, comparator: str = "less_than_or_equal") -> dict:
    return {"status": "pass", "residual": float(residual),
            "tolerance": float(tolerance), "comparator": comparator}


def fixture() -> tuple[dict, dict, dict]:
    src = gridding.gaussian_grid(8, 16, name="synthetic-atmosphere")
    dst = gridding.goldstein_grid(12, 8, igrid=gridding.GOLDSTEIN_EQUAL_AREA)
    src_contract = support_contract(src, "source_grid", "atmosphere_gaussian_grid")
    dst_contract = support_contract(dst, "destination_grid", "ocean_cgenie_grid")
    crossing = remap.Crossing(src, dst).restrict()
    rng = np.random.default_rng(20260828)
    probe = rng.normal(size=src.shape)
    extensive = probe * src.cell_area_fraction()
    energy = crossing.ledger(extensive, remap.EXTENSIVE_TOTAL)

    ones_dst, _ = crossing.apply(np.ones(src.shape), remap.INTENSIVE)
    same, _ = remap.Crossing(src, src).restrict().apply(probe, remap.INTENSIVE)
    back, _ = remap.Crossing(dst, src).restrict().apply(ones_dst, remap.INTENSIVE)
    reduced, _ = crossing.apply(probe, remap.INTENSIVE)
    processed_then, _ = crossing.apply(probe ** 2, remap.INTENSIVE)
    order_gap = float(np.max(np.abs(processed_then - reduced ** 2)))
    common_a, _ = remap.Crossing(src, src).restrict().apply(
        np.ones(src.shape), remap.INTENSIVE)
    common_b, _ = remap.Crossing(dst, src).restrict().apply(
        ones_dst, remap.INTENSIVE)

    na = {"status": "not_applicable", "reason": "synthetic fixture carries no such reservoir"}
    assessment = {
        "assessment_version": "vesper-spatial-conversion/1",
        "source_contract_identity": validate_contract(src_contract),
        "destination_contract_identity": validate_contract(dst_contract),
        "source_support_id": "source_grid",
        "destination_support_id": "destination_grid",
        "source_shape": list(src.shape),
        "destination_shape": list(dst.shape),
        "source_coordinates_sha256": src_contract["supports"][0]["coordinates_sha256"],
        "destination_coordinates_sha256": dst_contract["supports"][0]["coordinates_sha256"],
        "operator_version": "lib.remap.Crossing/1",
        "closure": {
            "area": pass_row(abs(src.cell_area_fraction().sum() - 1.0), 1e-14),
            "ocean_volume": copy.deepcopy(na),
            "water": copy.deepcopy(na),
            "salt": copy.deepcopy(na),
            "energy": pass_row(energy["closure_residual_relative"],
                               remap.CLOSURE_TOLERANCE),
            "carbon": copy.deepcopy(na),
            "nitrogen": copy.deepcopy(na),
            "phosphorus": copy.deepcopy(na),
        },
        "tests": {
            "constant": pass_row(float(np.max(np.abs(ones_dst - 1.0))),
                                 remap.CONSTANT_TOLERANCE),
            "identity": pass_row(float(np.max(np.abs(same - probe))), 1e-14),
            "reduction": pass_row(energy["closure_residual_relative"],
                                  remap.CLOSURE_TOLERANCE),
            "round_trip": pass_row(float(np.max(np.abs(back - 1.0))), 1e-14),
            "operator_order": pass_row(order_gap, 1e-6, "greater_than_or_equal"),
            "common_support": pass_row(float(np.max(np.abs(common_a - common_b))), 1e-14),
        },
        "inventory": {
            "nearest_fallbacks": [],
            "ownership_changes": [],
            "connectivity_changes": [],
            "discarded_spectral_content": [],
            "unmapped_extensive_stores": [],
        },
        "provenance": {
            "source_artifact_sha256": canonical_digest(probe.tolist()),
            "destination_geometry_sha256": dst_contract["supports"][0]["geometry_sha256"],
            "operator_sha256": canonical_digest({"operator": "lib.remap.Crossing/1"}),
        },
    }
    return assessment, src_contract, dst_contract


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=REPORT)
    args = parser.parse_args()
    spec = conversion_declaration()
    base, src, dst = fixture()
    checks: list[dict] = []

    def positive(name: str, candidate: dict) -> None:
        try:
            identity = validate_conversion_assessment(candidate, src, dst, spec)
            checks.append({"check": name, "pass": True, "identity": identity})
        except SpatialContractError as exc:
            checks.append({"check": name, "pass": False, "detail": str(exc)})

    def negative(name: str, mutate) -> None:
        candidate = copy.deepcopy(base)
        mutate(candidate)
        try:
            validate_conversion_assessment(candidate, src, dst, spec)
        except SpatialContractError as exc:
            checks.append({"check": name, "pass": True, "refusal": str(exc)})
        else:
            checks.append({"check": name, "pass": False,
                           "detail": "malformed assessment was accepted"})

    positive("complete assessed conversion", base)
    negative("source identity mismatch refuses",
             lambda c: c.update(source_contract_identity="0" * 64))
    negative("shape mismatch refuses", lambda c: c.update(source_shape=[8, 15]))
    negative("coordinate mismatch refuses",
             lambda c: c.update(destination_coordinates_sha256="0" * 64))
    negative("missing closure ledger refuses", lambda c: c["closure"].pop("salt"))
    negative("failed conserved quantity refuses",
             lambda c: c["closure"]["energy"].update(residual=1.0))
    negative("inapplicable closure needs a reason",
             lambda c: c["closure"]["water"].pop("reason"))
    negative("missing round-trip test refuses", lambda c: c["tests"].pop("round_trip"))
    negative("failed common-support test refuses",
             lambda c: c["tests"]["common_support"].update(residual=1.0))
    negative("operator-order test is discriminating",
             lambda c: c["tests"]["operator_order"].update(residual=0.0))
    negative("fallback inventory cannot disappear",
             lambda c: c["inventory"].pop("nearest_fallbacks"))
    negative("inventories stay machine-readable lists",
             lambda c: c["inventory"].update(ownership_changes="none"))
    negative("operator provenance is mandatory",
             lambda c: c["provenance"].pop("operator_sha256"))

    report = {
        "schema_version": spec["schema_version"],
        "assessment_version": spec["assessment_version"],
        "checks": checks,
        "fixture_assessment": base,
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
