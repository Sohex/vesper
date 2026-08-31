#!/usr/bin/env python3
"""Enforce SPAT-5's post-measurement representation decision."""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "exoplasim/scripts"))
CONTRACT = ROOT / "config/partial_surface.yaml"
REPORT = ROOT / "analysis/partial_surface_decision_gate_report.json"


def staged_field_populations(rung: str, nlat: int, nlon: int):
    """The three cell sets a tile's initial state is judged on, by INDEX.

    Both arrays are staged `.sra` on the same Gaussian grid and are read from
    the same directory, so the cells are matched by index and no longitude is
    reconstructed on either side. Code 172 is the binary ownership mask the
    model advances one surface on; code 1720 is the land fraction the tile
    model advances two on. Where they disagree is exactly the population whose
    initial state nothing has had to define until now.
    """
    from sra import read_sra
    staged = ROOT / "exoplasim/inputs" / rung.lower()
    fraction = read_sra(staged / f"orogen_{rung}_surf_1720.sra", nlat, nlon)
    owned = read_sra(staged / f"orogen_{rung}_surf_0172.sra", nlat, nlon)
    return staged, (owned < 0.5) & (fraction > 0.0), fraction >= 1.0


def uniform_over(field, population) -> bool:
    """Is the field a single value over this population?

    A field integrated from the native mesh varies over coastal cells on every
    continent, so one distinct value there is not a measurement of anything: it
    is the ownership mask's placeholder, written because the builder asked the
    binary mask where there was ground. That is the state the land tile would
    otherwise begin advancing from.
    """
    sub = field[population]
    return sub.size > 0 and np.unique(sub).size == 1


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

    # EACH MODULE BINARISES ITS OWN COPY, so there are three of these and not
    # one. landmod rounds `dls`, oceanmod `yls` and icemod `xls`, all at 0.5 and
    # all right after the mask is read. Checking oceanmod alone would report the
    # boundary as fractional while the land tile still saw a rounded mask, which
    # is the reading that costs the most: the seam would look taken.
    binarised = sorted(
        f"{unit}:{name}"
        for unit, name in (("landmod.f90", "dls"), ("oceanmod.f90", "yls"),
                           ("icemod.f90", "xls"))
        if all(fragment in (model_source / unit).read_text(
                   encoding="utf-8").lower()
               for fragment in (f"where ({name}(:) > 0.5)", f"{name}(:) = 1.0",
                                f"{name}(:) = 0.0")))
    hard_binary = bool(binarised)
    implementation_sources = {
        name: (model_source / name).read_text(encoding="utf-8").lower()
        for name in ("fluxmod.f90", "radmod.f90", "landmod.f90",
                     "seamod.f90", "outmod.f90")
    }
    # THE SEAMS ARE ORDERED, AND THE ORDER IS NOT A PREFERENCE. Exactly one
    # surface owner advances on a cell -- landstep is masked on dls > 0.0 and
    # seastep on dls < 0.5 -- and neither initialises outside its own mask. So
    # until both tiles carry their own state, the absent tile's boundary is a
    # COPY of the present one. An exchange seam taken first would evaluate the
    # same surface twice and area-weight two copies of it: a model that
    # reports per-tile fluxes which are not tiles, and which would then be
    # credited as implemented. STATE FIRST, EXCHANGE AFTER.
    # And the state seams themselves stand on a declaration. A tile that begins
    # advancing on a cell it did not own has to start from a value someone
    # named: without that, the seam invents an initial condition at the moment
    # it is switched on, which is the least inspectable place to put one.
    DECLARED_STATE = "declared initial state for both tiles"
    initial_state = (declaration.get("tile_model", {}).get("state", {})
                     .get("initial_state", {}))
    initial_state_fields = initial_state.get("fields", {})
    declared_available = sorted(
        name for name, entry in initial_state_fields.items()
        if entry.get("status") == "available")
    # AND `available` IS A CLAIM ABOUT THE ARTIFACT, NOT A STRING. Every field
    # here is one edit away from reading `available` while the builder that
    # would derive it has not run, and the gate would then credit the state
    # seam over a placeholder. So each field names the surface codes it reaches
    # the model through and the staged `.sra` is read: a field whose value is
    # uniform over the sea-owned cells that carry ground was decided by the
    # ownership mask and not by the mesh. An empty code list is a field the
    # model computes itself -- a cold-start profile or a scalar -- and there is
    # nothing staged to read.
    placeholder_fields: dict[str, list[int]] = {}
    verified_fields: dict[str, list[int]] = {}
    unstaged_fields: dict[str, list[int]] = {}
    planet = yaml.safe_load(
        (ROOT / "config/planet.yaml").read_text(encoding="utf-8"))["model"]
    rung = str(planet["resolution"]).upper()
    nlat, nlon = int(planet["latitudes"]), int(planet["longitudes"])
    staged_dir, sea_owned_with_ground, _pure_land = staged_field_populations(
        rung, nlat, nlon)
    from sra import read_sra
    for name in declared_available:
        codes = initial_state_fields[name].get("staged_surface_codes") or []
        for code in codes:
            path = staged_dir / f"orogen_{rung}_surf_{int(code):04d}.sra"
            if not path.is_file():
                unstaged_fields.setdefault(name, []).append(int(code))
                continue
            field = read_sra(path, nlat, nlon)
            target = (placeholder_fields if uniform_over(field, sea_owned_with_ground)
                      else verified_fields)
            target.setdefault(name, []).append(int(code))
    unready_initial_state = sorted(
        set(name for name, entry in initial_state_fields.items()
            if entry.get("status") != "available")
        | set(placeholder_fields) | set(unstaged_fields))
    check("every initial-state field declared available reaches the staged artifact",
          not placeholder_fields and not unstaged_fields,
          (f"{int(sea_owned_with_ground.sum())} sea-owned {rung} cells carry ground; "
           + (", ".join(f"{n} code {c} derived there"
                        for n, c in sorted(verified_fields.items()))
              or "no field declares a staged code")
           + (
               "; UNIFORM there, so the ownership mask decided it: "
               + ", ".join(f"{n} code {c}" for n, c in sorted(placeholder_fields.items()))
               if placeholder_fields else "")
           + (
               "; declared available but not staged: "
               + ", ".join(f"{n} code {c}" for n, c in sorted(unstaged_fields.items()))
               if unstaged_fields else "")))
    LAND_STATE = "separate positive-fraction land state"
    OCEAN_STATE = "separate positive-fraction ocean state"
    BOTH_STATES = (LAND_STATE, OCEAN_STATE)
    seam_order = (
        (DECLARED_STATE, None, (), ()),
        (LAND_STATE, "spat5_tile_state", ("landmod.f90",), (DECLARED_STATE,)),
        (OCEAN_STATE, "spat5_tile_state", ("seamod.f90",), (DECLARED_STATE,)),
        ("tile restart records", "spat5_tile_restart",
         ("landmod.f90", "seamod.f90"), BOTH_STATES),
        ("area-weighted turbulent exchange", "spat5_tile_turbulent",
         ("fluxmod.f90",), BOTH_STATES),
        ("area-weighted radiative exchange", "spat5_tile_radiative",
         ("radmod.f90",), BOTH_STATES),
        ("separate tile diagnostics", "spat5_tile_diagnostics",
         ("outmod.f90",), BOTH_STATES),
    )
    required_seams = {
        name: (bool(initial_state_fields) and not unready_initial_state
               if token is None
               else all(token in implementation_sources[unit] for unit in units))
        for name, token, units, _ in seam_order
    }
    missing_seams = [name for name, _, _, _ in seam_order
                     if not required_seams[name]]
    # A seam standing on an absent prerequisite is worse than a missing one:
    # it looks like progress and cannot be right.
    out_of_order = [
        name for name, _, _, prerequisites in seam_order
        if required_seams[name]
        and not all(required_seams[p] for p in prerequisites)
    ]
    next_seam = next(
        (name for name, _, _, prerequisites in seam_order
         if not required_seams[name]
         and all(required_seams[p] for p in prerequisites)),
        None)
    check("implementation seams are taken in dependency order",
          not out_of_order,
          ("no seam stands on an absent prerequisite" if not out_of_order
           else "state-dependent seams present without tile state: "
                + ", ".join(out_of_order)))
    check("selected representation is implemented",
          (not expected_tile or (not hard_binary and not missing_seams)),
          ("tile selected; "
           + (f"still hard-binarised at {', '.join(binarised)}; "
              if hard_binary else "")
           + (f"next seam is {next_seam}; missing " + ", ".join(missing_seams)
              if missing_seams else "all implementation seams present")))

    failed = [item for item in checks if not item["passed"]]
    report = {
        "contract_version": "vesper-partial-surface-decision-gate/1",
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "declaration": str(CONTRACT.relative_to(ROOT)),
        "measurement": str(measured_path.relative_to(ROOT)),
        "implementation": {
            "hard_binary_boundary": hard_binary,
            "hard_binarised_at": binarised,
            "required_seams": required_seams,
            "seam_order": [name for name, _, _, _ in seam_order],
            "initial_state_contract": initial_state.get("contract_version"),
            "initial_state_not_available": unready_initial_state,
            "initial_state_declared_available": declared_available,
            "sea_owned_cells_carrying_ground": int(sea_owned_with_ground.sum()),
            "staged_codes_derived_on_them": verified_fields,
            "staged_codes_uniform_on_them": placeholder_fields,
            "staged_codes_declared_but_absent": unstaged_fields,
            "seam_prerequisites": {name: list(prerequisites)
                                   for name, _, _, prerequisites in seam_order},
            "missing_seams": missing_seams,
            "out_of_order_seams": out_of_order,
            "next_seam": next_seam,
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
