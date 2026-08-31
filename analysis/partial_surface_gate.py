#!/usr/bin/env python3
"""SPAT-5 declaration, geometry report and flux-bracket regression gate."""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
sys.path.insert(0, str(ROOT / "lib"))

import provenance                                             # noqa: E402
import rungs                                                  # noqa: E402
from coastline_flux_bracket import bracket_field              # noqa: E402

REPORT = ROOT / "analysis/partial_surface_gate_report.json"
# The surface code the fractional land support is carried on. SPAT-5.
LAND_FRACTION_CODE = 1720


def tile_combine(fraction: np.ndarray, land: np.ndarray,
                 ocean: np.ndarray) -> np.ndarray:
    """Endpoint-preserving extensive combination used by SPAT-5 fixtures.

    Zero-area tiles are not evaluated: multiplying an uninitialised/poisoned
    absent tile by zero would still propagate NaN and violate the declared
    bit-identical pure-surface reductions.
    """
    fraction = np.asarray(fraction, dtype=np.float64)
    land = np.asarray(land, dtype=np.float64)
    ocean = np.asarray(ocean, dtype=np.float64)
    if fraction.shape != land.shape or land.shape != ocean.shape:
        raise ValueError("fraction, land and ocean arrays must have one shape")
    if np.any((fraction < 0.0) | (fraction > 1.0)):
        raise ValueError("tile fraction must be in [0, 1]")
    out = np.empty_like(fraction)
    pure_land = fraction == 1.0
    pure_ocean = fraction == 0.0
    partial = ~(pure_land | pure_ocean)
    out[pure_land] = land[pure_land]
    out[pure_ocean] = ocean[pure_ocean]
    equal_tiles = partial & (land == ocean)
    mixed_tiles = partial & ~equal_tiles
    out[equal_tiles] = land[equal_tiles]
    out[mixed_tiles] = (fraction[mixed_tiles] * land[mixed_tiles]
                        + (1.0 - fraction[mixed_tiles]) * ocean[mixed_tiles])
    return out


def _routine_body(source: str, name: str) -> str:
    """One Fortran routine's body, from lowercased source.

    Splitting on the bare name is wrong and quietly returns the WRONG region:
    `end subroutine <name>` contains `subroutine <name>`, so the last part is
    everything AFTER the routine and every substring test on it passes or
    fails for reasons that have nothing to do with the routine. The opening
    line is anchored to its own newline and indentation instead.
    """
    opening = f"\n      subroutine {name}\n"
    if opening not in source:
        return ""
    body = source.split(opening, 1)[1]
    return body.split(f"\n      end subroutine {name}", 1)[0]


def main() -> None:
    checks: list[dict] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})

    cfg = yaml.safe_load((ROOT / "config/planet.yaml").read_text())
    contract = yaml.safe_load((ROOT / "config/partial_surface.yaml").read_text())
    area = json.loads((ROOT / "analysis/coastline_threshold_cost.json").read_text())
    # The rung is read from the configuration, never spelled: the staged inputs
    # are keyed by it and a gate that names one rung passes on a tree running
    # another while reading a directory that is not the one the model starts
    # from. SPAT-2.
    rung, _, _ = rungs.model_grid(cfg)
    inputs = ROOT / "exoplasim" / "inputs" / rung.lower()
    boundary = json.loads((inputs / "boundary_conditions_report.json").read_text())
    check("active build identity", area["source_build"] == cfg["source_build"],
          f"{area['source_build']} vs {cfg['source_build']}")
    check("whole ladder measured",
          [r["rung"] for r in area["rungs"]] == ["T21", "T42", "T85", "T127", "T170"],
          ", ".join(r["rung"] for r in area["rungs"]))
    check("class partitions close", all(
        r["ledger"]["class_partition_max_residual"] <= 1e-12 for r in area["rungs"]),
        "maximum residual " + str(max(
            r["ledger"]["class_partition_max_residual"] for r in area["rungs"])))
    check("third class is absent, not small",
          area["third_class"]["mesh_area_fraction"] == 0.0
          and area["third_class"]["status"] == "absent_not_small",
          area["third_class"]["status"])
    check("selection state is exclusive and measured",
          contract["selection"]["strategy"] == "cost_bracket_first_decide_after"
          and bool(contract["selection"]["tile_model_selected"])
              != bool(contract["selection"]["model_form_bracket_selected"]),
          contract["selection"]["tile_model_status"])
    # `staged_surface_field` and not a constructed path: it is the one door onto
    # a staged `.sra`, it refuses a field staged from another build, and it
    # refuses an unstamped one. A gate REPORTS rather than aborts, so its
    # refusal becomes this check's detail.
    try:
        staged = provenance.staged_surface_field(LAND_FRACTION_CODE, cfg)
        staged_detail = f"staged from {staged['build']}"
    except SystemExit as refusal:
        staged, staged_detail = None, str(refusal)
    # The boundary report is written for the configured rung, so the ledger row
    # it must agree with is that rung's, chosen by name and not by position.
    measured = next((r for r in area["rungs"] if r["rung"] == rung), None)
    check("fractional land support is emitted",
          LAND_FRACTION_CODE in boundary.get("codes", [])
          and boundary.get("land_fraction", {}).get("surface_code") == LAND_FRACTION_CODE
          and measured is not None
          and boundary.get("land_fraction", {}).get("partial_cells") ==
              measured["ledger"]["partial_land_cells"]
          and staged is not None,
          (f"code {boundary.get('land_fraction', {}).get('surface_code')}, "
           f"{boundary.get('land_fraction', {}).get('partial_cells')} partial "
           f"{rung} cells; {staged_detail}"))

    model_state = (ROOT / "vendor/exoplasim/exoplasim/plasim/src/"
                   "plasimmod.f90").read_text(encoding="utf-8").lower()
    check("compiled endpoint-preserving tile combine exists",
          "elemental real function spat5_tile_combine" in model_state
          and "if (pfraction == 1.0)" in model_state
          and "else if (pfraction == 0.0)" in model_state,
          "shared primitive for turbulent, radiative and diagnostic routing")

    model_driver = (ROOT / "vendor/exoplasim/exoplasim/plasim/src/"
                    "plasim.f90").read_text(encoding="utf-8").lower()
    surface_source = (ROOT / "vendor/exoplasim/exoplasim/plasim/src/"
                      "surfmod.f90").read_text(encoding="utf-8").lower()
    tile_boundary_records = (
        "dlt_ts", "dlt_qs", "dlt_rhs", "dlt_z0", "dlt_alb", "dlt_sa1",
        "dlt_sa2", "dot_ts", "dot_qs", "dot_rhs", "dot_z0", "dot_alb",
        "dot_sa1", "dot_sa2")
    check("tile boundary archives are restart-complete",
          all(f"call mpputgp('{name}'" in model_driver
              and f"call mpgetgp('{name}'" in model_driver
              for name in tile_boundary_records)
          and "call spat5_archive_land_boundary" in surface_source
          and "call spat5_archive_ocean_boundary" in surface_source
          and "call spat5_complete_absent_tile" in surface_source,
          f"{len(tile_boundary_records)} land/ocean boundary records")

    # The absent tile must be completed over the EXECUTION mask. Keyed on the
    # fraction's endpoints instead, the completion covers only pure cells and
    # leaves every partial cell holding landmod's 1.0e20 sentinel or seamod's
    # unassigned array -- which is the whole population the tile model is for,
    # and it reaches the restart. Measured at T21: 620 and 433 of 1053.
    absent_tile_body = _routine_body(surface_source,
                                     "spat5_complete_absent_tile")
    check("absent tile is completed over the execution mask",
          "where (dls(:) > 0.5)" in absent_tile_body
          and "elsewhere" in absent_tile_body
          and "dlf(:) == 0.0" not in absent_tile_body
          and "dlf(:) == 1.0" not in absent_tile_body,
          "keyed on dls, so partial cells are covered and no sentinel survives")

    # It must run before anything can read the bundles on a continuation:
    # fluxstep and radstep precede the first surfstep, so a restart-only
    # completion leaves a legacy sentinel live for one timestep.
    surfini_body = _routine_body(surface_source, "surfini")
    check("absent-tile completion is unconditional at initialisation",
          "call spat5_complete_absent_tile" in surfini_body
          and "if (nrestart == 0)    call spat5_complete_absent_tile"
              not in surfini_body,
          "cold start and continuation both repair the absent tile")
    flux_source = (ROOT / "vendor/exoplasim/exoplasim/plasim/src/"
                   "fluxmod.f90").read_text(encoding="utf-8").lower()
    radiation_source = (ROOT / "vendor/exoplasim/exoplasim/plasim/src/"
                        "radmod.f90").read_text(encoding="utf-8").lower()
    check("surface owners consume separate seeded exchange channels",
          "call spat5_load_land_exchange" in surface_source
          and "call spat5_load_ocean_exchange" in surface_source
          and "call spat5_restore_exchange_aggregate" in surface_source
          and "call spat5_seed_turbulent_bundles" in flux_source
          and "dlt_swfl(:)=dswfl(:,nlep)" in radiation_source
          and "dot_swfl(:)=dswfl(:,nlep)" in radiation_source,
          "land/ocean turbulent+radiative bundles; dual evaluation pending")

    # Endpoint branches are part of the operator, not an optimisation. Poison
    # the absent tile so a multiply-by-zero implementation fails visibly.
    tf = np.array([1.0, 0.0, 0.25, 0.75])
    tl = np.array([17.0, np.nan, 4.0, 8.0])
    to = np.array([np.nan, 23.0, 12.0, 0.0])
    combined = tile_combine(tf, tl, to)
    check("tile pure-surface reductions bypass absent state",
          np.array_equal(combined[:2], np.array([17.0, 23.0])),
          str(combined[:2]))
    check("tile extensive combination closes",
          np.array_equal(combined[2:], np.array([10.0, 6.0])),
          str(combined[2:]))
    signed_zero = tile_combine(np.array([0.37]), np.array([-0.0]),
                               np.array([-0.0]))
    check("identical tile values preserve their bit pattern",
          signed_zero.view(np.uint64)[0] == np.array([-0.0]).view(np.uint64)[0],
          f"0x{signed_zero.view(np.uint64)[0]:016x}")
    try:
        tile_combine(np.array([-np.finfo(float).eps]), np.array([1.0]),
                     np.array([1.0]))
    except ValueError:
        bounds_refused = True
    else:
        bounds_refused = False
    check("tile fraction bounds fail closed", bounds_refused,
          "negative one-epsilon fixture refused" if bounds_refused else "accepted")

    # Two rows, with four cells of each class globally. Partial cells at [0,0]
    # and [1,3] need the opposite class. The first uses same-row land values
    # 10 and 14; the second uses same-row ocean values 2 and 4.
    field = np.array([[2., 4., 10., 14.], [2., 4., 10., 14.]])
    land = np.array([[False, False, True, True],
                     [False, False, True, True]])
    fraction = np.array([[0.25, 0., 1., 1.],
                         [0., 0., 1., 0.75]])
    cell_area = np.ones_like(field)
    row = bracket_field(field, fraction, land, cell_area, 2)
    check("synthetic lower bracket",
          np.isclose(row["global_mean_delta_lower"], -1.0 / 8.0),
          str(row["global_mean_delta_lower"]))
    check("synthetic upper bracket",
          np.isclose(row["global_mean_delta_upper"], 0.5 / 8.0),
          str(row["global_mean_delta_upper"]))
    check("pure cells contribute exactly zero", row["partial_cells"] == 2,
          f"{row['partial_cells']} partial cells")
    check("extensive closure", row["extensive_closure_residual"] <= 1e-14,
          str(row["extensive_closure_residual"]))

    # With one target cell in the starting row and minimum two, the reference
    # search must expand and inventory that fact rather than silently switching
    # populations.
    expanded = bracket_field(field, fraction, land, cell_area, 3)
    check("row expansion is explicit",
          expanded["reference_inventory"]["expanded_cells"] == 2
          and expanded["reference_inventory"]["maximum_row_expansion"] == 1,
          json.dumps(expanded["reference_inventory"], sort_keys=True))

    failed = [x for x in checks if not x["passed"]]
    report = {
        "contract_version": "vesper-partial-surface-gate/1",
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "checks": checks,
        "summary": {"passed": len(checks) - len(failed), "failed": len(failed)},
        "verdict": "PASS" if not failed else "FAIL",
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n")
    for item in checks:
        print(f"{'PASS' if item['passed'] else 'FAIL'}  {item['name']}: {item['detail']}")
    print(f"\n{report['verdict']}: {len(checks)-len(failed)} passed, {len(failed)} failed; "
          f"{REPORT.relative_to(ROOT)}")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
