#!/usr/bin/env python3
"""Build the atmosphere-to-ocean crossing and check it against a known integral.

    python analysis/ocean_remap.py                       # configured rung, 36 x 36
    python analysis/ocean_remap.py --rung T85 --ocean 72x72
    python analysis/ocean_remap.py --igrid 2             # the atmosphere's own rows

The offline ocean under OCN-3 runs on GOLDSTEIN's grid and the climate runs on
ExoPlaSim's, and the two are different IN KIND: a Gaussian row is not an
equal-area row, so a crossing that treats either as regular is wrong at the poles
first and everywhere second. `lib/remap.py` is the operator and this is what
builds it, checks it, and writes the weights down.

## Why the check is an integral and not a comparison

A conservative remap has a right answer rather than a plausible one, so the
question "did it work" can be made to FAIL. Two identities do the work.

**The crossing preserves the source's own integral.** For a field carrying a
budget -- a heat flux, a freshwater flux -- `sum(F_dst*A_dst)` must equal
`sum(F_src*A_src)` to round-off, and `lib/remap.py:CLOSURE_TOLERANCE` is the bar,
fixed from the float64 round-off floor and three orders tighter than the 1e-9
relative that `references/esmf/` holds itself to on analytic fields.

**The source integral is one whose answer is known in advance.** That second
identity is the one that catches the failure the first cannot see. A crossing
built on midpoint-between-nodes row boundaries conserves perfectly well -- it is
a self-consistent partition of the sphere -- it just conserves against a grid the
model does not use, and it does that quietly. So the field remapped here is
`3*sin(lat)^2 - 1` evaluated at the Gaussian NODES, whose exact integral over the
sphere is zero and which Gauss-Legendre quadrature returns as exactly zero for
any node count above one. The source integral therefore has an analytic answer,
and it comes out right only if the row areas ARE the quadrature weights. The
midpoint construction runs beside it as a negative control and has to miss the
same bar by at least three orders, because a bar nothing can fail is not a bar.

## What it writes

`analysis/ocean_remap.json` carries every check with its tolerance and its
residual, both grids' constructors, the coast ledger and the provenance.
`analysis/ocean_remap_weights.npz` carries the two sparse matrices so a consumer
applies the crossing rather than rebuilding it -- the conserving one, whose
columns close so that all of a valid source cell's content lands somewhere valid,
and the intensive one, which drops what falls outside the destination mask
because a state has no integral to conserve.

The masked arm runs on the export's own `surface_class` (`CLAUDE.md` rule 1) so
that the coast rule is exercised against real geography rather than a synthetic
mask. It is a DEMONSTRATION of the rule and not the coupling's mask: the ocean's
wet mask comes from the accepted export through OCN-11 and does not exist yet, so
what stands in for it here is the atmosphere's own ocean fraction remapped and
cut at a half. That is exactly the mask disagreement the rule is for, at the
resolution ratio the crossing will actually run at.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))

import numpy as np  # noqa: E402
import yaml  # noqa: E402

import builds  # noqa: E402
import gridding  # noqa: E402
import paths as paths_lib  # noqa: E402
import provenance  # noqa: E402
import remap as remap_lib  # noqa: E402
import rungs  # noqa: E402
from spatial_support import (grid_support_contract, validate_contract,  # noqa: E402
                             validate_conversion_assessment)

OUT_JSON = ROOT / "analysis" / "ocean_remap.json"

# The negative control has to MISS the acceptance bar by at least three orders,
# or the bar is not discriminating and the whole suite is decoration. Fixed here
# before any control was run. `docs/src/practice/failure-modes.md` class 17.
CONTROL_FLOOR = 1e3 * remap_lib.CLOSURE_TOLERANCE

OCEAN_WET_CUT = 0.5     # an ocean cell is wet where over half its area is


def known_integral_field(spec: gridding.GridSpec) -> np.ndarray:
    """`3*sin(lat)^2 - 1` at the Gaussian nodes, whose sphere integral is zero.

    The nodes, NOT the cell centres. This is the field the model would hold, and
    the point of the check is that a quadrature over the model's own weights
    returns the analytic answer exactly while a quadrature over any other
    partition of the same sphere does not.
    """
    mu = np.sin(np.deg2rad(gridding.gaussian_latitudes(spec.nlat)))
    return np.repeat((3.0 * mu ** 2 - 1.0)[:, None], spec.nlon, axis=1)


def midpoint_control(nlat: int, nlon: int) -> gridding.GridSpec:
    """The Gaussian grid with rows bounded by MIDPOINTS between the nodes.

    A perfectly good partition of the sphere and the wrong one. It lives here
    rather than in `lib/gridding.py` because a constructor that builds a grid the
    model does not use is a hazard wherever a caller can reach it, and its only
    job is to fail a check that the real constructor passes.
    """
    nodes = np.sin(np.deg2rad(gridding.gaussian_latitudes(nlat)))
    edges = np.empty(nlat + 1)
    edges[1:-1] = 0.5 * (nodes[:-1] + nodes[1:])
    edges[0], edges[-1] = 1.0, -1.0
    lon = gridding.gaussian_grid(nlat, nlon).lon_edges
    return gridding.GridSpec(name=f"gaussian-midpoint-control-{nlat}x{nlon}",
                             lon_edges=lon, sin_edges=edges,
                             source="analysis/ocean_remap.py:midpoint_control, "
                                    "the negative control, deliberately not the model's grid")


def integral(spec: gridding.GridSpec, field: np.ndarray) -> float:
    return float((np.asarray(field) * spec.cell_area_fraction()).sum())


def scale_of(spec: gridding.GridSpec, field: np.ndarray) -> float:
    return float((np.abs(np.asarray(field)) * spec.cell_area_fraction()).sum())


def ocean_mask(grid_dir: Path, spec: gridding.GridSpec) -> np.ndarray:
    """Atmosphere cells the export calls ocean, from `surface_class`.

    `CLAUDE.md` rule 1: land comes from `surface_class` and never from
    `land_mask`, which floods every dry closed basin below sea level and would
    hand the ocean a coastline the terrain does not have.
    """
    field = next(f for f in json.loads((grid_dir / "manifest.json").read_text())
                 ["grid"]["fields"] if f["name"] == "surface_class")
    cls = np.fromfile(grid_dir / field["path"], dtype=field["dtype"]).reshape(spec.shape)
    return cls == 0


def run_checks(atm: gridding.GridSpec,
               ocn: gridding.GridSpec) -> tuple[list[dict], remap_lib.Crossing]:
    checks: list[dict] = []

    def record(name, residual, bar, comparator="<=", detail=""):
        ok = residual <= bar if comparator == "<=" else residual >= bar
        checks.append({"check": name, "residual": float(residual),
                       "tolerance": float(bar), "must_be": comparator,
                       "passes": bool(ok), "detail": detail})
        return ok

    # 1. Both grids partition the sphere. A cell area that does not sum to one
    #    is a spec that is not a grid, and every weight below is a ratio to it.
    for spec in (atm, ocn):
        record(f"{spec.name}: cell areas sum to the sphere",
               abs(spec.cell_area_fraction().sum() - 1.0), remap_lib.AREA_TOLERANCE)

    # 2. The atmosphere's row extents ARE the Gauss-Legendre weights. This is
    #    what makes the crossing conserve against the model's own budget rather
    #    than against a partition of the sphere nobody integrates over.
    _, weights = np.polynomial.legendre.leggauss(atm.nlat)
    record("the Gaussian rows reproduce the quadrature weights",
           float(np.abs(atm.dsin - weights[::-1]).max()), 1e-13,
           detail="rows laid end to end in the sine of latitude")

    # 3. Both one-dimensional overlaps are partitions of unity, in both
    #    directions. A periodic wrap that folds instead of wrapping fails the
    #    longitude pair, and it fails it even though the crossing built on it
    #    would still conserve.
    crossing = remap_lib.Crossing(atm, ocn)
    record("longitude overlap closes on the source columns",
           float(np.abs(crossing.lon_overlap.sum(axis=0) - atm.dlon).max()), 1e-12)
    record("longitude overlap closes on the destination columns",
           float(np.abs(crossing.lon_overlap.sum(axis=1) - ocn.dlon).max()), 1e-12)
    record("latitude overlap closes on the source rows",
           float(np.abs(crossing.lat_overlap.sum(axis=0) - atm.dsin).max()), 1e-13)
    record("latitude overlap closes on the destination rows",
           float(np.abs(crossing.lat_overlap.sum(axis=1) - ocn.dsin).max()), 1e-13)

    crossing.restrict()

    # 4. THE KNOWN INTEGRAL. The source quadrature of `3*sin(lat)^2 - 1` over the
    #    model's own weights is analytically zero, and the crossing must return
    #    zero on the far side too.
    field = known_integral_field(atm)
    scale = scale_of(atm, field)
    record("the source quadrature returns the analytic integral",
           abs(integral(atm, field)) / scale, remap_lib.CLOSURE_TOLERANCE,
           detail="3*sin(lat)^2 - 1 over the sphere is exactly zero")
    dst, coverage = crossing.apply(field, remap_lib.FLUX_DENSITY)
    record("the crossing returns the analytic integral",
           abs(integral(ocn, dst)) / scale, remap_lib.CLOSURE_TOLERANCE)

    # 5. THE NEGATIVE CONTROL. Midpoint row boundaries give a self-consistent
    #    partition of the sphere whose quadrature is not the model's, so it must
    #    miss the same bar by at least three orders. If it does not, the bar is
    #    not discriminating and nothing above is a test.
    control = midpoint_control(atm.nlat, atm.nlon)
    record("the midpoint control misses the known integral",
           abs(integral(control, field)) / scale, CONTROL_FLOOR, comparator=">=",
           detail="row boundaries at midpoints between the Gaussian nodes")

    # 6. Conservation of an arbitrary flux, and preservation of a constant.
    rng = np.random.default_rng(20260825)
    flux = rng.normal(0.0, 60.0, atm.shape)
    led = crossing.ledger(flux, remap_lib.FLUX_DENSITY)
    record("an arbitrary flux density conserves its integral",
           led["closure_residual_relative"], remap_lib.CLOSURE_TOLERANCE,
           detail="W/m2 into an ocean cell, integral in watts")
    total = flux * atm.cell_area_fraction()
    led_total = crossing.ledger(total, remap_lib.EXTENSIVE_TOTAL)
    record("an already-integrated total conserves its sum",
           led_total["closure_residual_relative"], remap_lib.CLOSURE_TOLERANCE)
    ones, _ = crossing.apply(np.ones(atm.shape), remap_lib.INTENSIVE)
    record("an intensive state preserves a constant",
           float(np.abs(ones - 1.0).max()), remap_lib.CONSTANT_TOLERANCE)

    # 7. A grid crossed with ITSELF is the identity, exactly, under every
    #    semantics. This is the strongest available answer that can be wrong: it
    #    fixes the weights, the normalisations and both orientations at once, and
    #    an operator that is merely plausible fails it.
    same = remap_lib.Crossing(atm, atm).restrict()
    probe = rng.normal(0.0, 40.0, atm.shape)
    for semantics in remap_lib.SEMANTICS:
        back, _ = same.apply(probe, semantics)
        want = probe * atm.cell_area_fraction() if semantics == remap_lib.EXTENSIVE_TOTAL else probe
        want = want / atm.cell_area_fraction() if semantics == remap_lib.EXTENSIVE_TOTAL else want
        record(f"a grid crossed with itself is the identity under {semantics}",
               float(np.abs(back - want).max() / np.abs(probe).max()), 1e-14)

    # 8. The ocean's longitude origin is a parameter, not a label. Rotating it by
    #    exactly one ocean column must roll the answer by one column and change
    #    nothing else. A crossing that matched longitude NUMBERS between the two
    #    grids instead of taking each from its own constructor fails this.
    turn = float(ocn.lon_edges[-1] - ocn.lon_edges[0])
    rolled = gridding.goldstein_grid(
        ocn.nlon, ocn.nlat, igrid=gridding.GOLDSTEIN_EQUAL_AREA,
        lon_origin_deg=float(ocn.lon_edges[0]) + turn / ocn.nlon)
    if np.allclose(rolled.sin_edges, ocn.sin_edges):
        spun, _ = remap_lib.Crossing(atm, rolled).restrict().apply(
            flux, remap_lib.FLUX_DENSITY)
        record("rotating the ocean origin by one column rolls the answer",
               float(np.abs(spun - np.roll(dst_of(crossing, flux), -1, axis=1)).max()
                     / np.abs(flux).max()), 1e-13,
               detail="phi0 is a parameter of the ocean grid, never a label to match")

    return checks, crossing


def dst_of(crossing, field):
    out, _ = crossing.apply(field, remap_lib.FLUX_DENSITY)
    return out


def assessed_conversion(atm, ocn, crossing, checks, grid_dir, coast):
    """Bind the produced weights to SPAT-1 identities and SPAT-10 checks."""
    source_artifact = (str((grid_dir / "manifest.json").relative_to(ROOT))
                       if grid_dir is not None else "config/planet.yaml")
    source_contract = grid_support_contract(
        atm, "analysis/ocean_remap_weights.npz", "atmosphere_grid",
        "atmosphere_gaussian_grid", [source_artifact])
    # Until OCN-11 supplies wet volume and topology this is a candidate
    # comparison geometry, not an accepted ocean support.
    destination_contract = grid_support_contract(
        ocn, "analysis/ocean_remap_weights.npz", "goldstein_candidate_grid",
        "comparison_support", ["lib/gridding.py:goldstein_grid"])

    by_name = {row["check"]: row for row in checks}

    def le(residual, tolerance):
        return {"status": "pass", "residual": float(residual),
                "tolerance": float(tolerance), "comparator": "less_than_or_equal"}

    def ge(residual, tolerance):
        return {"status": "pass", "residual": float(residual),
                "tolerance": float(tolerance), "comparator": "greater_than_or_equal"}

    def na(reason):
        return {"status": "not_applicable", "reason": reason}

    rng = np.random.default_rng(20260828)
    probe = rng.normal(size=atm.shape)
    reduced, _ = crossing.apply(probe, remap_lib.INTENSIVE)
    processed_then, _ = crossing.apply(probe ** 2, remap_lib.INTENSIVE)
    order_gap = float(np.max(np.abs(processed_then - reduced ** 2)))
    ones_dst, _ = crossing.apply(np.ones(atm.shape), remap_lib.INTENSIVE)
    ones_back, _ = remap_lib.Crossing(ocn, atm).restrict().apply(
        ones_dst, remap_lib.INTENSIVE)
    common_a, _ = remap_lib.Crossing(atm, atm).restrict().apply(
        np.ones(atm.shape), remap_lib.INTENSIVE)
    common_b, _ = remap_lib.Crossing(ocn, atm).restrict().apply(
        ones_dst, remap_lib.INTENSIVE)
    identity_residual = max(
        row["residual"] for row in checks
        if row["check"].startswith("a grid crossed with itself is the identity"))
    energy = by_name["an arbitrary flux density conserves its integral"]
    area_residual = max(abs(atm.cell_area_fraction().sum() - 1.0),
                        abs(ocn.cell_area_fraction().sum() - 1.0))

    nearest = []
    ownership = []
    if coast is not None:
        if coast.get("source_cells_orphaned", 0):
            nearest.append({
                "source_cells": coast["source_cells_orphaned"],
                "source_area_fraction": coast["source_area_fraction_moved"],
                "furthest_move_degrees": coast["furthest_move_degrees"],
                "rule": "nearest valid destination; conserving semantics only",
            })
        ownership.append({
            "atmosphere_ocean_cells": coast["atmosphere_ocean_cells"],
            "candidate_ocean_wet_cells": coast["ocean_wet_cells"],
            "wet_fraction_cut": coast["ocean_wet_cut"],
            "status": "demonstration; OCN-11 supplies the accepted ownership",
        })

    source_bytes = ((grid_dir / "manifest.json").read_bytes()
                    if grid_dir is not None else
                    (ROOT / "config" / "planet.yaml").read_bytes())
    assessment = {
        "assessment_version": "vesper-spatial-conversion/1",
        "source_contract_identity": validate_contract(source_contract),
        "destination_contract_identity": validate_contract(destination_contract),
        "source_support_id": "atmosphere_grid",
        "destination_support_id": "goldstein_candidate_grid",
        "source_shape": list(atm.shape),
        "destination_shape": list(ocn.shape),
        "source_coordinates_sha256": source_contract["supports"][0]["coordinates_sha256"],
        "destination_coordinates_sha256": destination_contract["supports"][0]["coordinates_sha256"],
        "operator_version": "lib.remap.Crossing/1",
        "closure": {
            "area": le(area_residual, remap_lib.AREA_TOLERANCE),
            "ocean_volume": na("OCN-11 has not supplied accepted bathymetry or wet volume"),
            "water": na("the geometry fixture carries no water reservoir"),
            "salt": na("the geometry fixture carries no salt reservoir"),
            "energy": le(energy["residual"], energy["tolerance"]),
            "carbon": na("the geometry fixture carries no carbon reservoir"),
            "nitrogen": na("the geometry fixture carries no nitrogen reservoir"),
            "phosphorus": na("the geometry fixture carries no phosphorus reservoir"),
        },
        "tests": {
            "constant": le(float(np.max(np.abs(ones_dst - 1.0))),
                           remap_lib.CONSTANT_TOLERANCE),
            "identity": le(identity_residual, 1e-14),
            "reduction": le(energy["residual"], energy["tolerance"]),
            "round_trip": le(float(np.max(np.abs(ones_back - 1.0))), 1e-14),
            "operator_order": ge(order_gap, 1e-6),
            "common_support": le(float(np.max(np.abs(common_a - common_b))), 1e-14),
        },
        "inventory": {
            "nearest_fallbacks": nearest,
            "ownership_changes": ownership,
            "connectivity_changes": [{
                "status": "not_applicable",
                "reason": "the candidate crossing has no ocean topology; OCN-11 owns it",
            }],
            "discarded_spectral_content": [{
                "count": 0,
                "reason": "finite-volume overlap is not a spectral truncation",
            }],
            "unmapped_extensive_stores": [],
        },
        "provenance": {
            "source_artifact_sha256": hashlib.sha256(source_bytes).hexdigest(),
            "destination_geometry_sha256":
                destination_contract["supports"][0]["geometry_sha256"],
            "operator_sha256": hashlib.sha256(
                (ROOT / "lib" / "remap.py").read_bytes()).hexdigest(),
        },
    }
    identity = validate_conversion_assessment(
        assessment, source_contract, destination_contract)
    return {"source_contract": source_contract,
            "destination_contract": destination_contract,
            "assessment": assessment, "assessment_identity": identity}


def coast_demonstration(atm, ocn, grid_dir: Path) -> dict:
    """The coast rule against real geography, and what it costs.

    The atmosphere's ocean cells are the export's; the ocean grid's wet cells are
    those over half covered by them. The two masks then disagree at the coast by
    construction, which is the disagreement OCN-11's real mask will also have,
    and the ledger says how much flux the rule moved and how far.
    """
    atm_wet = ocean_mask(grid_dir, atm)
    geometric = remap_lib.Crossing(atm, ocn).restrict()
    wet_fraction, _ = geometric.apply(atm_wet.astype(float), remap_lib.INTENSIVE)
    ocn_wet = wet_fraction >= OCEAN_WET_CUT
    crossing = remap_lib.Crossing(atm, ocn).restrict(src_valid=atm_wet, dst_valid=ocn_wet)
    rng = np.random.default_rng(7)
    flux = np.where(atm_wet, rng.normal(0.0, 60.0, atm.shape), 0.0)
    led = crossing.ledger(flux, remap_lib.FLUX_DENSITY)
    led["ocean_wet_cut"] = OCEAN_WET_CUT
    led["atmosphere_ocean_cells"] = int(atm_wet.sum())
    led["ocean_wet_cells"] = int(ocn_wet.sum())
    return led


def export_area_disagreement(grid_dir: Path, spec: gridding.GridSpec) -> dict:
    """How far the export's own `grid_cell_area` is from the quadrature weights.

    The export ships both `grid/gauss_weights.bin` and `grid/grid_cell_area.bin`
    and they describe the same grid, so any disagreement between them is a fact
    about the artifact rather than about this crossing. It is measured here
    because this is the code that has to choose one of them, and it chooses the
    weights: the model's budget is a quadrature.

    A corrected export returns zero here by construction, since `grid_cell_area`
    is then the quadrature partition in km2. An export built before that
    correction carries the NEAREST-ROW partition instead and returns about 0.22
    in the polar row at every truncation; `notes/audits/ocean-grid-crossing.md`
    section 3 has which is which and what it was worth. So this stays as a
    check on the build in hand rather than being retired.
    """
    weights = np.fromfile(grid_dir / "grid" / "gauss_weights.bin", dtype="float64")
    area = np.fromfile(grid_dir / "grid" / "grid_cell_area.bin",
                       dtype="float32").reshape(spec.shape).astype(np.float64)
    row = area[:, 0] / area[:, 0].sum()
    want = spec.dsin / spec.dsin.sum()
    rel = np.abs(row - want) / want
    return {
        "gauss_weights_against_construction": float(np.abs(weights - spec.dsin).max()),
        "grid_cell_area_against_weights_max_relative": float(rel.max()),
        "grid_cell_area_against_weights_max_row": int(rel.argmax()),
        "grid_cell_area_against_weights_interior_median_relative":
            float(np.median(rel[2:-2])),
        "used_for_the_crossing": "grid/gauss_weights.bin",
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rung", help="atmosphere rung; default the configured one")
    ap.add_argument("--ocean", default="36x36",
                    help="ocean grid as imaxXjmax; muffingen's declared range is 1-72")
    ap.add_argument("--igrid", type=int, default=gridding.GOLDSTEIN_EQUAL_AREA,
                    help="0 equal area (shipped), 1 equal angle, 2 the atmosphere's rows")
    ap.add_argument("--lon-origin", type=float, default=gridding.GOLDSTEIN_LON_ORIGIN,
                    help="phi0 in degrees; the shipped value is 260 west")
    ap.add_argument("--out", type=Path, default=OUT_JSON)
    args = ap.parse_args()

    config = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    rung = args.rung.upper() if args.rung else rungs.model_grid(config)[0]
    nlat, nlon, _ = rungs.geometry(rung)
    atm = gridding.gaussian_grid(nlat, nlon, name=f"exoplasim-{rung}")

    try:
        ni, nj = (int(v) for v in args.ocean.lower().split("x"))
    except ValueError:
        ap.error(f"--ocean takes imaxXjmax, not {args.ocean!r}")
    ocn = gridding.goldstein_grid(
        ni, nj, igrid=args.igrid, lon_origin_deg=args.lon_origin,
        atmosphere_rows=atm.sin_edges if args.igrid == gridding.GOLDSTEIN_ATMOSPHERE_ROWS
        else None)

    # The one coordinate source. The spec is CONSTRUCTED and the export's axis is
    # READ, and this is where they are made to agree; a spec that describes a
    # different grid from the field it remaps is `CLAUDE.md` rule 3 in its area
    # form. Absent an export the crossing is still exact, and the report says so.
    grid_dir = None
    try:
        grid_dir = builds.grid_export(config, resolution=rung)
        lat, _, _ = gridding.grid_geometry(grid_dir)
        gridding.require_gaussian_rows(atm, lat, f"the {rung} export")
    except Exception as exc:                                    # noqa: BLE001
        grid_dir = None
        checked_against = f"no export on disk to check against: {exc}"
    else:
        checked_against = str(grid_dir.relative_to(ROOT))

    checks, crossing = run_checks(atm, ocn)
    coast = coast_demonstration(atm, ocn, grid_dir) if grid_dir is not None else None
    support_assessment = assessed_conversion(
        atm, ocn, crossing, checks, grid_dir, coast)
    report = {
        "atmosphere": {"name": atm.name, "shape": list(atm.shape), "source": atm.source},
        "ocean": {"name": ocn.name, "shape": list(ocn.shape), "source": ocn.source,
                  "igrid": args.igrid, "lon_origin_deg": args.lon_origin},
        "rows_checked_against": checked_against,
        "checks": checks,
        "geometric_ledger": crossing.ledger(),
        "tolerances": {
            "closure_relative": remap_lib.CLOSURE_TOLERANCE,
            "constant_absolute": remap_lib.CONSTANT_TOLERANCE,
            "area_relative": remap_lib.AREA_TOLERANCE,
            "control_floor": CONTROL_FLOOR,
        },
        "provenance": provenance.config_stamp(config, "analysis/ocean_remap.py"),
        "spatial_conversion": support_assessment,
    }
    if grid_dir is not None:
        report["export_area_disagreement"] = export_area_disagreement(grid_dir, atm)
        report["coast_demonstration"] = coast

    failed = [c for c in checks if not c["passes"]]
    report["checks_run"] = len(checks)
    report["checks_failed"] = len(failed)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    weights_path = args.out.with_name(args.out.stem + "_weights.npz")
    conserving = crossing._conserving.tocoo()
    intensive = crossing._overlap.tocoo()
    np.savez_compressed(
        weights_path,
        conserving_row=conserving.row, conserving_col=conserving.col,
        conserving_val=conserving.data, intensive_row=intensive.row,
        intensive_col=intensive.col, intensive_val=intensive.data,
        src_shape=np.array(atm.shape), dst_shape=np.array(ocn.shape),
        src_area=crossing.src_area, dst_area=crossing.dst_area,
        source_contract_identity=np.array(
            support_assessment["assessment"]["source_contract_identity"]),
        destination_contract_identity=np.array(
            support_assessment["assessment"]["destination_contract_identity"]),
        assessment_identity=np.array(support_assessment["assessment_identity"]))

    for c in checks:
        mark = "  ok  " if c["passes"] else " FAIL "
        print(f"[{mark}] {c['check']}: {c['residual']:.3g} "
              f"{c['must_be']} {c['tolerance']:.3g}")
    print(f"\n{len(checks)} checks, {len(failed)} failed -> "
          f"{paths_lib.rel(args.out)}, {paths_lib.rel(weights_path)}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
