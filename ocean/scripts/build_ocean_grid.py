#!/usr/bin/env python3
"""Place the mesh on the ocean's grid, cross a vector onto it, and check both.

    python ocean/scripts/build_ocean_grid.py
    python ocean/scripts/build_ocean_grid.py --ocean 72x72 --rung T85
    python ocean/scripts/build_ocean_grid.py --igrid 1

The offline ocean under OCN-3 runs on GOLDSTEIN's grid, which no exporter
writes: it exists as a constructor in `lib/gridding.py` and nowhere else. Two
crossings onto it are therefore not the ones the tree already has, and this
builds and checks both.

**The mesh onto the ocean grid.** `lib/gridding.py:region_cells` bins the mesh
onto a grid the export SHIPS and takes the row boundaries from that grid's own
axis. There is no such axis here, so a mesh region is placed against the
`GridSpec`'s cell BOUNDARIES by `spec_cells`, in the sine of latitude because
the ocean's rows are uniform in the sine and not in the angle. This is what
OCN-11's wet area, volume and topology will be reduced over, and it is the
reason the ocean's own bathymetry has to be written through this same placement:
the frame the ocean model runs in is then the generator's by construction rather
than by agreement.

**A vector onto the ocean grid, and back off it.** The ocean returns a surface
VELOCITY to the next ExoPlaSim baseline, which world-pt8 advects sea ice with.
"East" at one longitude is not "east" at another, so remapping the two
components as independent scalars averages numbers that are components in
different frames. `lib/remap.py:Crossing.apply_vector` lifts the pair into the
sphere's own three cartesian components at the source cells' frames, remaps
those with the same weights, and projects onto the destination cells' frames.

## The two checks, and what would mean wrong

**The vector crossing.** A rigid rotation `omega x r` has an analytically zero
cartesian integral over the sphere, because the integral is `omega` crossed into
the integral of position and that is the origin. Its exact cell mean is
`omega x c` at the cell's own area centroid, so the source integral is zero
before any remapping and a conserving crossing has to return zero on the far
side.

THAT INTEGRAL IS NOT BY ITSELF A TEST, and the control is what says so. Both
grids are uniform in longitude, so the frame error is a pure phase and cancels
around every row: the COMPONENTWISE remap returns the analytic zero to
round-off as well. This runs that case deliberately and asserts the control
passes it, because a reader otherwise takes the full-sphere integral for a
verdict on the lift. The discriminating form is the same conservation law over
a region that is NOT zonally symmetric -- here a 90 degree longitude wedge
placed on the rotation axis, so nothing about it is chosen against a grid --
where the source's own integral is still `omega` crossed into a closed-form sum
of cell centroids, the vector crossing reproduces it, and the componentwise
remap misses it by four orders. The real ocean mask is not zonally symmetric
either, which is why this is the case that matters rather than the tidy one.

**The placement.** Every mesh region placed by `spec_cells` lies between the
boundaries of the cell it landed in, taken from the GOLDSTEIN constructor's own
edge arrays rather than from a restatement of them, and the area ledger closes.
Beside it, the same door run on the atmosphere's Gaussian grid is compared with
`region_cells`, which is the door for that grid: the two disagree by
construction, because a Gaussian node is a quadrature abscissa and not the
centre of its cell, and the disagreement is measured here so that nobody reaches
for the wrong door on the strength of it looking like the right one.

## What it writes

`ocean/data/<build>/ocean_grid.json` carries the ocean grid's SPAT-1 support
contract and identity, every check with its tolerance and residual, the
placement ledger and resolution diagnostics, and the provenance.
`ocean/data/<build>/ocean_grid_placement.npz` carries the per-region flat cell
index, so OCN-11 reduces over the placement rather than rebuilding it.

## What this is NOT

It is not the ocean support. Wet area and volume, partial coasts, hypsometry,
shelves and sills, routed river mouths and bathymetric smoothing are OCN-11's
and are not decided here; `ocean/scripts/build_spatial_support.py` owns them.
Nothing in this file declares a cell wet. The GOLDSTEIN geometry stays a
candidate comparison support, exactly as `analysis/ocean_remap.py` labels it,
until that contract exists.

It writes no NetCDF, and that is a decision rather than an omission.
`lib/nc_geometry.py` declares a product's geometry in one of the project's two
longitude conventions, both of which are namings of the Gaussian grid the
atmosphere and the export share. The ocean grid's origin is `phi0`, a rotation
of a different grid, and a file declaring it would need a third arm in that
module. The consumer that needs one is OCN-11's `spatial_support.nc`, so the arm
belongs in the commit that writes it and not ahead of it.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import yaml

from _paths import DATA, PLANET_CONFIG, PROJECT_ROOT

import builds                                                    # noqa: E402
import gridding                                                  # noqa: E402
import orogen                                                    # noqa: E402
import paths as paths_lib                                        # noqa: E402
import provenance                                                # noqa: E402
import remap as remap_lib                                        # noqa: E402
import rungs                                                     # noqa: E402
from spatial_support import grid_support_contract, validate_contract  # noqa: E402

# The rotation axis, in the equatorial plane and at a longitude that is not a
# column boundary of either grid, so no check passes on a symmetry it was handed.
ROTATION_AXIS_LON_DEG = 37.0

# The discriminating check is taken over a longitude window rather than the whole
# sphere, because the whole sphere is where the frame error cancels. A quarter
# turn, placed on the axis above so the window is chosen against the field and
# not against either grid's columns. A destination cell joins it when this much
# of its area is covered, which is `analysis/ocean_remap.py`'s wet cut and the
# same argument: a cell more than half covered has valid source under it.
WEDGE_WIDTH_DEG = 90.0
WEDGE_COVERAGE_CUT = 0.5

# Below this share of the peak mapped speed a per-cell relative error is a ratio
# of two small numbers and says nothing, so the diagnostic states its population
# rather than reporting a quotient that grows without bound near the field's own
# null. `docs/src/practice/failure-modes.md` class 34.
FIELD_FLOOR_FRACTION = 0.2

# THE ACCEPTANCE BARS, FIXED BEFORE ANY RESULT WAS SEEN.
#
# The vector integral is a sum over cells of a quantity that cancels to zero, so
# its floor is `lib/remap.py`'s own closure floor and it takes the same bar,
# whose derivation from float64 round-off is stated there.
VECTOR_TOLERANCE = remap_lib.CLOSURE_TOLERANCE

# The control has to MISS that bar by at least three orders or the bar is not
# discriminating and the check above is decoration. Same floor and same argument
# as `analysis/ocean_remap.py:CONTROL_FLOOR`.
CONTROL_FLOOR = 1e3 * remap_lib.CLOSURE_TOLERANCE

# The placement is exact set membership, so its only tolerance is the round-off
# of taking a sine and folding a longitude by a whole turn. A degree of
# longitude at float64 carries a few ulp of 360; 1e-9 is orders above that and
# orders below the narrowest column either grid has.
PLACEMENT_SLACK = 1e-9

# The bar on the mesh's own longitude reproducing from its cartesian
# coordinates. The export stores both in float32, whose ulp at 180 degrees is
# about 1e-05, and the round trip goes through an arctangent; 1e-03 is two
# orders above that and four below the narrowest column either grid has, which
# is the disagreement it would have to catch.
COORDINATE_SLACK_DEG = 1e-3


def _igrid_of(spec: gridding.GridSpec) -> int:
    """The `igrid` arm a GOLDSTEIN spec was built at, off its own constructor.

    `GridSpec.source` records what constructed the edges, so the arm is read
    back from the spec rather than passed alongside it and allowed to disagree.
    """
    for arm in (gridding.GOLDSTEIN_EQUAL_AREA, gridding.GOLDSTEIN_EQUAL_ANGLE,
                gridding.GOLDSTEIN_ATMOSPHERE_ROWS):
        if f"igrid={arm}" in spec.source:
            return arm
    raise ValueError(f"{spec.name} was not built by goldstein_grid: {spec.source}")


def check_row(name, residual, bar, comparator="<=", detail=""):
    ok = residual <= bar if comparator == "<=" else residual >= bar
    return {"check": name, "residual": float(residual), "tolerance": float(bar),
            "must_be": comparator, "passes": bool(ok), "detail": detail}


def rotation_field(spec: gridding.GridSpec, axis: np.ndarray) -> np.ndarray:
    """`omega x r` as the EXACT cell mean over each cell, shaped (ncell, 3).

    A rigid rotation is linear in position, so its mean over a cell is the
    rotation applied to the cell's mean position, which `cell_centroids` gives
    analytically. Taking it at a cell centre instead would be an approximation
    and the source integral would then be zero only to the order of a cell.

    The field is perpendicular to the centroid and so lies in the tangent plane
    of the frame built on it, which is why the east/north pair below carries all
    of it and the source radial component is zero rather than small.
    """
    return np.cross(axis, spec.cell_centroids())


def components(spec: gridding.GridSpec, cart: np.ndarray):
    """A cartesian field as (east, north) on the grid, plus what did not fit."""
    east, north, up = spec.cell_frames()
    return ((cart * east).sum(axis=1).reshape(spec.shape),
            (cart * north).sum(axis=1).reshape(spec.shape),
            float(np.abs((cart * up).sum(axis=1)).max()))


def integral(spec: gridding.GridSpec, cart: np.ndarray) -> np.ndarray:
    """The area-weighted cartesian integral of a per-cell 3-vector field."""
    return (cart * spec.cell_area_fraction().ravel()[:, None]).sum(axis=0)


def wedge(spec: gridding.GridSpec, start_deg: float, width_deg: float) -> np.ndarray:
    """Source cells inside a longitude window, shaped like the grid.

    A window and not the whole sphere, because the whole sphere is where the
    frame error cancels. Its start is the ROTATION AXIS's own longitude, so the
    window is placed against the field rather than against either grid's columns
    and no check is handed a boundary that happens to line up.
    """
    _, centres = spec.cell_centres()
    lon = np.tile(centres, spec.nlat)
    folded = lon - 360.0 * np.floor((lon - start_deg) / 360.0)
    return (folded < start_deg + width_deg).reshape(spec.shape)


def vector_checks(atm: gridding.GridSpec, ocn: gridding.GridSpec) -> tuple[list[dict], dict]:
    """The rigid rotation, the crossing, and the control that has to miss."""
    e_d, n_d, u_d = ocn.cell_frames()
    lam = np.deg2rad(ROTATION_AXIS_LON_DEG)
    axis = np.array([np.cos(lam), np.sin(lam), 0.0])
    src = rotation_field(atm, axis)
    east, north, src_radial = components(atm, src)
    src_area = atm.cell_area_fraction().ravel()
    speed = np.linalg.norm(src, axis=1)
    scale = float((speed * src_area).sum())

    def mapped(crossing, e, n):
        E, N, _, R = crossing.apply_vector(e, n, remap_lib.FLUX_DENSITY)
        return (E.ravel()[:, None] * e_d + N.ravel()[:, None] * n_d
                + R.ravel()[:, None] * u_d), R

    def componentwise(crossing, e, n):
        ce, _ = crossing.apply(e, remap_lib.FLUX_DENSITY)
        cn, _ = crossing.apply(n, remap_lib.FLUX_DENSITY)
        return ce.ravel()[:, None] * e_d + cn.ravel()[:, None] * n_d

    rows = [
        check_row("the rigid rotation is tangent at every source cell's own frame",
                  src_radial / scale, VECTOR_TOLERANCE,
                  detail="omega x c is perpendicular to c, so the east/north pair loses nothing"),
        check_row("the source quadrature returns the analytic integral",
                  float(np.abs(integral(atm, src)).max()) / scale, VECTOR_TOLERANCE,
                  detail="the integral of omega x r over the sphere is omega crossed into the origin"),
    ]

    whole = remap_lib.Crossing(atm, ocn).restrict()
    vec, radial = mapped(whole, east, north)
    rows.append(check_row(
        "the vector crossing returns the analytic integral over the sphere",
        float(np.abs(integral(ocn, vec)).max()) / scale, VECTOR_TOLERANCE,
        detail="lifted to cartesian at the source frames, remapped, projected at the destination's"))
    comp = componentwise(whole, east, north)
    rows.append(check_row(
        "the full-sphere integral does NOT separate the crossing from the control",
        float(np.abs(integral(ocn, comp)).max()) / scale, VECTOR_TOLERANCE,
        detail="both grids are uniform in longitude, so the frame error is a phase that "
               "cancels around every row; asserted so the zero above is not read as a verdict"))

    # THE DISCRIMINATING FORM: the same conservation law over a region that is
    # not zonally symmetric. The source's own integral over the window is still
    # closed form -- omega crossed into a sum of analytic cell centroids -- and
    # `restrict` closes the columns, so all of a valid source cell's content
    # lands somewhere valid and the identity survives the mask.
    src_valid = wedge(atm, ROTATION_AXIS_LON_DEG, WEDGE_WIDTH_DEG)
    coverage, _ = whole.apply(src_valid.astype(float), remap_lib.INTENSIVE)
    dst_valid = coverage >= WEDGE_COVERAGE_CUT
    masked = remap_lib.Crossing(atm, ocn).restrict(src_valid=src_valid, dst_valid=dst_valid)
    want = np.cross(axis, (atm.cell_centroids()[src_valid.ravel()]
                           * src_area[src_valid.ravel(), None]).sum(axis=0))
    wedge_scale = float((speed * src_area)[src_valid.ravel()].sum())
    e_w = np.where(src_valid, east, 0.0)
    n_w = np.where(src_valid, north, 0.0)
    vec_w, _ = mapped(masked, e_w, n_w)
    rows.append(check_row(
        "the vector crossing returns the analytic integral over a wedge",
        float(np.abs(integral(ocn, vec_w) - want).max()) / wedge_scale, VECTOR_TOLERANCE,
        detail=f"a {WEDGE_WIDTH_DEG:.0f} degree window on the rotation axis, where the "
               "frame error has nothing to cancel against"))
    comp_w = componentwise(masked, e_w, n_w)
    rows.append(check_row(
        "the componentwise control misses the wedge integral",
        float(np.abs(integral(ocn, comp_w) - want).max()) / wedge_scale,
        CONTROL_FLOOR, comparator=">=",
        detail="east and north remapped as two independent scalars, which averages "
               "components in frames that are not the same frame"))

    rows.append(check_row(
        "the crossing discards no vector where the two grids agree",
        float(np.abs(remap_lib.Crossing(ocn, ocn).restrict().apply_vector(
            *components(ocn, rotation_field(ocn, axis))[:2],
            remap_lib.FLUX_DENSITY)[3]).max()) / scale,
        VECTOR_TOLERANCE,
        detail="a grid crossed with itself has one tangent plane, so the radial part is zero"))

    # WHAT THE LIFT IS WORTH, in the units world-pt8 cares about: how far the
    # componentwise answer is from the crossing's, per cell, against the local
    # speed rather than a global scale. Reported below a fifth of the peak speed
    # is a ratio of two small numbers, so the population is cut there and the
    # cut is stated.
    got = np.linalg.norm(vec, axis=1)
    gap = np.linalg.norm(vec - comp, axis=1)
    live = got > FIELD_FLOOR_FRACTION * got.max()
    diagnostics = {
        "frame_error_over_local_speed": {
            "population": "destination cells above "
                          f"{FIELD_FLOOR_FRACTION:.2f} of the peak mapped speed",
            "cells": int(live.sum()),
            "median": float(np.median(gap[live] / got[live])),
            "p99": float(np.quantile(gap[live] / got[live], 0.99)),
            "max": float((gap[live] / got[live]).max()),
        },
        "radial_discarded_over_local_speed_max":
            float(np.abs(radial).ravel()[live].max() / got[live].min()),
        "wedge": {
            "start_deg": ROTATION_AXIS_LON_DEG, "width_deg": WEDGE_WIDTH_DEG,
            "coverage_cut": WEDGE_COVERAGE_CUT,
            "source_cells": int(src_valid.sum()),
            "destination_cells": int(dst_valid.sum()),
            "orphans": masked.orphans,
        },
    }
    return rows, diagnostics


def placement(ocn: gridding.GridSpec, export, radius_km: float) -> dict:
    """Place the mesh on the ocean grid and say what the placement carries."""
    lat = np.asarray(export.lat, dtype=np.float64)
    lon = np.asarray(export.lon, dtype=np.float64)
    area = np.asarray(export.cell_area, dtype=np.float64)
    cell = gridding.spec_cells(ocn, lat, lon)

    row, col = np.divmod(cell, ocn.nlon)
    sine = np.sin(np.deg2rad(lat))
    turn = float(ocn.lon_edges[-1] - ocn.lon_edges[0])
    folded = lon - turn * np.floor((lon - ocn.lon_edges[0]) / turn)
    lo_s = np.minimum(ocn.sin_edges[row], ocn.sin_edges[row + 1])
    hi_s = np.maximum(ocn.sin_edges[row], ocn.sin_edges[row + 1])
    outside = int(((sine < lo_s - PLACEMENT_SLACK)
                   | (sine > hi_s + PLACEMENT_SLACK)
                   | (folded < ocn.lon_edges[col] - PLACEMENT_SLACK)
                   | (folded > ocn.lon_edges[col + 1] + PLACEMENT_SLACK)).sum())

    # RULE 3, CHECKED RATHER THAN INHERITED. The placement compares a mesh
    # coordinate against edges built from the ocean model's own setup, and that
    # is sound only if the two are in ONE frame. Three things say so, and each
    # of them can fail.
    #
    # First, `lon` is the generator's own coordinate and not a label: it
    # reproduces from the mesh's cartesian `x, y, z`, which are y-up, so the
    # polar axis is `y` and the longitude is `atan2(x, z)`. A label read off a
    # file reproduces from nothing.
    #
    # Second, the fold into the spec's window is a WHOLE number of turns, which
    # is a rotation of nothing: the same physical longitude, renumbered. It is
    # the same shift `lib/remap.py:_periodic_overlap` applies between the same
    # two constructors, and a reconciliation that moved a region by any other
    # amount would show up here as a residue.
    #
    # Third, `phi0` is a parameter of the ocean grid and not a naming of the
    # atmosphere's: rotating it by exactly one ocean column has to move every
    # region exactly one column and change nothing else. A placement that read
    # one grid's column number as the other's cannot do that, and this is
    # `analysis/ocean_remap.py`'s origin check in the mesh's form.
    reproduced = np.rad2deg(np.arctan2(np.asarray(export.x, dtype=np.float64),
                                       np.asarray(export.z, dtype=np.float64)))
    coordinate_gap = float(np.abs(gridding.longitude_difference(lon, reproduced)).max())
    fold_residue = float(np.abs(gridding.longitude_difference(folded, lon)).max())
    rolled = gridding.goldstein_grid(
        ocn.nlon, ocn.nlat, igrid=_igrid_of(ocn),
        lon_origin_deg=float(ocn.lon_edges[0]) + turn / ocn.nlon,
        atmosphere_rows=ocn.sin_edges if _igrid_of(ocn) ==
        gridding.GOLDSTEIN_ATMOSPHERE_ROWS else None)
    rolled_col = gridding.spec_cells(rolled, lat, lon) % ocn.nlon
    origin_moved = int((rolled_col != (col - 1) % ocn.nlon).sum())

    ledger = gridding.transfer_ledger(cell, ocn.ncell, area)
    counts = np.bincount(cell, minlength=ocn.ncell)
    landed = np.bincount(cell, weights=area, minlength=ocn.ncell)
    want = ocn.cell_area_fraction().ravel() * (4.0 * np.pi * radius_km ** 2)
    ratio = landed / want
    return {
        "cell": cell,
        "checks": [
            check_row("every mesh region lies between the boundaries of the cell it landed in",
                      outside, 0.0,
                      detail=f"{export.n_regions} regions against the GOLDSTEIN "
                             "constructor's own sine and longitude edges"),
            check_row("the placement loses no mesh area",
                      ledger["closure_residual_relative"], remap_lib.CLOSURE_TOLERANCE,
                      detail="the mesh region dual areas offered against those that landed"),
            check_row("no ocean cell is empty of mesh",
                      int((counts == 0).sum()), 0.0,
                      detail="an empty cell is a cell whose content would have to be invented"),
            check_row("the mesh longitude reproduces from the mesh's own cartesian coordinates",
                      coordinate_gap, COORDINATE_SLACK_DEG,
                      detail="atan2(x, z) on a y-up frame; a label read off a file "
                             "reproduces from nothing, which is what makes this a frame "
                             "and not a naming"),
            check_row("the fold into the ocean's window is a whole number of turns",
                      fold_residue, PLACEMENT_SLACK,
                      detail="the same shift lib/remap.py:_periodic_overlap applies, and a "
                             "rotation of nothing rather than a reconciliation"),
            check_row("rotating phi0 by one ocean column moves every region one column",
                      origin_moved, 0.0,
                      detail="phi0 is a parameter of the ocean grid, never a label to match "
                             "against the atmosphere's"),
            check_row("the mesh area in a cell is the cell's own area",
                      float(np.abs(ratio - 1.0).max()), 0.02,
                      detail="the mesh is a discrete sample of the sphere, so this closes to "
                             "the sampling noise of the regions in the smallest cell"),
        ],
        "ledger": ledger,
        "regions_per_cell": {
            "min": int(counts.min()), "median": float(np.median(counts)),
            "max": int(counts.max()), "mean": float(counts.mean()),
        },
        "mesh_area_over_cell_area": {
            "min": float(ratio.min()), "max": float(ratio.max()),
            "median": float(np.median(ratio)),
        },
    }


def door_disagreement(export, grid_dir: Path, spec: gridding.GridSpec) -> dict:
    """What the WRONG door costs on the atmosphere's own grid, measured.

    `region_cells` is the door for a grid the export ships and bins to the row
    the export itself used, which is the nearest Gaussian NODE. `spec_cells`
    bins between the quadrature EDGES. A node is not the centre of its cell, so
    the two put a band of mesh either side of every row boundary in different
    rows. The columns are the same expression in both and must agree exactly; if
    they ever did not, the disagreement below would be rule 3 rather than this.
    """
    by_axis, _, nlon = gridding.region_cells(export, grid_dir)
    by_edges = gridding.spec_cells(spec, np.asarray(export.lat, dtype=np.float64),
                                   np.asarray(export.lon, dtype=np.float64))
    area = np.asarray(export.cell_area, dtype=np.float64)
    moved = by_axis != by_edges
    same_column = int((by_axis % nlon != by_edges % nlon).sum())
    return {
        "grid": spec.name,
        "regions_placed_differently": int(moved.sum()),
        "mesh_area_fraction_placed_differently": float(area[moved].sum() / area.sum()),
        "columns_placed_differently": same_column,
        "reading": "region_cells is the door for a grid the export ships; this is what "
                   "reaching for the constructed-boundary door there would move",
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rung", help="atmosphere rung for the vector crossing; "
                                   "default the configured one")
    ap.add_argument("--ocean", default="36x36",
                    help="ocean grid as imaxXjmax; muffingen's declared range is 1-72")
    ap.add_argument("--igrid", type=int, default=gridding.GOLDSTEIN_EQUAL_AREA,
                    help="0 equal area (shipped), 1 equal angle, 2 the atmosphere's rows")
    ap.add_argument("--lon-origin", type=float, default=gridding.GOLDSTEIN_LON_ORIGIN,
                    help="phi0 in degrees; the shipped value is 260 west")
    ap.add_argument("--no-mesh", action="store_true",
                    help="grid identity and the vector crossing only, with no export read")
    ap.add_argument("--out", type=Path,
                    help="report path; default ocean/data/<build>/ocean_grid.json")
    args = ap.parse_args()

    config = yaml.safe_load(PLANET_CONFIG.read_text(encoding="utf-8"))
    rung = args.rung.upper() if args.rung else rungs.model_grid(config)[0]
    nlat, nlon, _ = rungs.geometry(rung)
    atm = gridding.gaussian_grid(nlat, nlon, name=f"exoplasim-{rung}")

    try:
        ni, nj = (int(v) for v in args.ocean.lower().split("x"))
    except ValueError:
        return ap.error(f"--ocean takes imaxXjmax, not {args.ocean!r}")
    ocn = gridding.goldstein_grid(
        ni, nj, igrid=args.igrid, lon_origin_deg=args.lon_origin,
        atmosphere_rows=atm.sin_edges if args.igrid == gridding.GOLDSTEIN_ATMOSPHERE_ROWS
        else None)

    build = provenance.active_build(config)
    out = args.out or (DATA / build / "ocean_grid.json")

    checks, vector_diagnostics = vector_checks(atm, ocn)
    report = {
        "ocean_grid": {"name": ocn.name, "shape": list(ocn.shape),
                       "source": ocn.source, "igrid": args.igrid,
                       "lon_origin_deg": args.lon_origin},
        "atmosphere_grid": {"name": atm.name, "shape": list(atm.shape),
                            "source": atm.source, "rung": rung},
        "rotation_axis_longitude_deg": ROTATION_AXIS_LON_DEG,
        "source_build": build,
        "vector_crossing": vector_diagnostics,
    }

    placed = None
    if args.no_mesh:
        report["mesh"] = {"status": "not read", "reason": "--no-mesh"}
    else:
        mesh_dir = builds.mesh_export(config)
        export = orogen.Export(mesh_dir)
        radius_km = export.radius_km
        placed = placement(ocn, export, radius_km)
        checks += placed["checks"]
        grid_dir = builds.grid_export(config, resolution=rung)
        lat_axis, _, _ = gridding.grid_geometry(grid_dir)
        gridding.require_gaussian_rows(atm, lat_axis, f"the {rung} export")
        report["mesh"] = {
            "export": str(mesh_dir.relative_to(PROJECT_ROOT)),
            "terrain_hash": export.terrain_hash,
            "regions": int(export.n_regions),
            "radius_km": radius_km,
            "rows_checked_against": str(grid_dir.relative_to(PROJECT_ROOT)),
            "ledger": placed["ledger"],
            "regions_per_ocean_cell": placed["regions_per_cell"],
            "mesh_area_over_cell_area": placed["mesh_area_over_cell_area"],
            "wrong_door_on_the_atmosphere_grid":
                door_disagreement(export, grid_dir, atm),
        }

    contract = grid_support_contract(
        ocn, str(out.relative_to(PROJECT_ROOT)), "goldstein_candidate_grid",
        "comparison_support", ["lib/gridding.py:goldstein_grid"])
    report["ocean_grid"]["support_contract"] = contract
    report["ocean_grid"]["support_contract_identity"] = validate_contract(contract)
    report["ocean_grid"]["acceptance"] = (
        "a candidate comparison support, not an accepted ocean support: wet mask, "
        "volume and topology are OCN-11's and no line here declares a cell wet")

    report["tolerances"] = {
        "vector_relative": VECTOR_TOLERANCE,
        "control_floor": CONTROL_FLOOR,
        "placement_slack_degrees": PLACEMENT_SLACK,
    }
    report["checks"] = checks
    report["checks_run"] = len(checks)
    failed = [c for c in checks if not c["passes"]]
    report["checks_failed"] = len(failed)
    report["provenance"] = provenance.config_stamp(
        config, "ocean/scripts/build_ocean_grid.py")

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    written = [out]
    if placed is not None:
        npz = out.with_name(out.stem + "_placement.npz")
        np.savez_compressed(
            npz, region_cell=placed["cell"].astype(np.int32),
            ocean_shape=np.array(ocn.shape),
            lon_edges=ocn.lon_edges, sin_edges=ocn.sin_edges,
            support_contract_identity=np.array(
                report["ocean_grid"]["support_contract_identity"]),
            terrain_hash=np.array(report["mesh"]["terrain_hash"]))
        written.append(npz)

    for c in checks:
        mark = "  ok  " if c["passes"] else " FAIL "
        print(f"[{mark}] {c['check']}: {c['residual']:.3g} "
              f"{c['must_be']} {c['tolerance']:.3g}")
    print(f"\n{len(checks)} checks, {len(failed)} failed -> "
          + ", ".join(paths_lib.rel(p) for p in written))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
