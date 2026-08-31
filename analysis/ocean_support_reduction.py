#!/usr/bin/env python3
"""Size the ocean support crossing's nonlinear reductions on terrain alone.

    python analysis/ocean_support_reduction.py                  # the ladder and 36x36
    python analysis/ocean_support_reduction.py --rungs T42
    python analysis/ocean_support_reduction.py --no-connectivity

SPAT-7's ocean half. `notes/audits/ocean-and-marine-biosphere.md` section 8b
says narrow connections and partial coastal cells can change circulation even
where the cell-mean depth is unchanged, and section 8d says the coastal
boundary's geometry is what any delivery law would read. Neither sizes
anything. This does, on terrain alone: bathymetry and coastline come out of the
Orogen export and need no climate run and no ocean run.

`notes/audits/ocean-support-nonlinear-reductions.md` carries the finding and
the pre-registration -- the bars, the invariants, and the statement that no sign
is assumed for refinement -- and it was committed before this script ran.

## The three orders being compared

A cell's bathymetry is a DISTRIBUTION over the mesh regions the cell holds, and
three different quantities can be asked of it.

**Volume** is affine in depth, so the area-weighted mean reproduces it exactly.
That arm is the CONTROL: it can only fail through a binning or an area error,
and `lib/remap.py:CLOSURE_TOLERANCE` is what it is judged against.

**The area on either side of a depth horizon** is a threshold on the
distribution, which the mean does not carry. `lib/gridding.py:cell_expectation`
is the operator -- evaluate the indicator per region and average it -- and the
gap against testing the cell mean is what a shelf area, a hypsometric class or
a deep-water partition is built on.

**Connectivity** is neither: it is a property of the graph the cells form, and
it moves in two directions at once. A cell mean is at least the cell minimum,
so coarsening DEEPENS bottlenecks and opens passages; a passage narrower than a
cell can fall below the ocean fraction the atmosphere binarises at, which closes
it. The measurement is which dominates and at what depth.

## What it does not do

No ocean model runs, and none is needed: every arm here is geometry. The
passage-capacity bar is therefore reported as the critical fractional error
against an ocean heat transport the first accepted ocean run supplies, rather
than collapsed onto a guessed circulation. `CLAUDE.md` rule 3 holds throughout:
the Gaussian rungs are binned with `gridding.region_cells`, which is the binning
every surface builder uses, and the GOLDSTEIN grid with a searchsorted binner on
the spec's own constructed edges. No longitude label is matched across any
boundary.
"""

from __future__ import annotations

import argparse
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
import rungs  # noqa: E402
from orogen import Export, OCEAN  # noqa: E402

OUT_JSON = ROOT / "analysis" / "ocean_support_reduction.json"
TRANSPORT_LOOP = ROOT / "ocean" / "config" / "transport_loop.yaml"
PARTIAL_SURFACE = ROOT / "config" / "partial_surface.yaml"

# The horizons the arms are evaluated at, in km below the datum. Declared here
# rather than derived from the bathymetry, because a horizon chosen from the
# distribution it is about to cut is a horizon chosen after the result.
HORIZONS_KM = tuple(round(0.25 * k, 3) for k in range(0, 25))

# The exactness bar for the volume control. `lib/remap.py` fixes it as the
# float64 round-off floor every conservative crossing in this project closes
# against, and this arm closes the same kind of integral.
VOLUME_TOLERANCE = 1.0e-12

SEA_WATER_DENSITY_KG_M3 = 1030.0    # only for the printed capacity conversion


def _criterion(loop: dict, name: str) -> dict:
    for c in loop["exit"]["criteria"]:
        if c["id"] == name:
            return c
    raise SystemExit(f"ocean/config/transport_loop.yaml has no criterion {name!r}")


def spec_cells(spec: gridding.GridSpec, lat_deg, lon_deg) -> np.ndarray:
    """Flat cell index on a constructed spec, for points in export coordinates.

    The spec's own edges are the only coordinate source: the row comes from the
    SINE of latitude against `sin_edges` in whichever order that array runs, and
    the column from `lon_edges` after wrapping into the spec's own origin. The
    check in `main` is that this reproduces `gridding.column` exactly on a grid
    whose columns are the export's, which is what says the wrap is the same
    convention rather than a second opinion about it.
    """
    s = np.sin(np.deg2rad(np.asarray(lat_deg, dtype=np.float64)))
    edges = spec.sin_edges
    if edges[0] < edges[-1]:
        r = np.searchsorted(edges, s, side="right") - 1
    else:
        r = np.searchsorted(-edges, -s, side="right") - 1
    r = np.clip(r, 0, spec.nlat - 1)
    origin = float(spec.lon_edges[0])
    wrapped = np.mod(np.asarray(lon_deg, dtype=np.float64) - origin, 360.0) + origin
    c = np.clip(np.searchsorted(spec.lon_edges, wrapped, side="right") - 1,
                0, spec.nlon - 1)
    return r.astype(np.int64) * spec.nlon + c.astype(np.int64)


def grid_edges(nlat: int, nlon: int) -> tuple[np.ndarray, np.ndarray]:
    """Undirected 4-neighbour edges on an (nlat, nlon) grid, periodic in longitude.

    A grid ocean model's cells exchange across their shared faces, so this is the
    graph connectivity is asked of. Longitude wraps and latitude does not, which
    is the sphere's own topology on a latitude-longitude mesh; the poles are not
    joined across, because no cell pair there shares a face.
    """
    rows = np.arange(nlat)[:, None]
    cols = np.arange(nlon)[None, :]
    flat = rows * nlon + cols
    east = np.stack([flat.ravel(), (rows * nlon + (cols + 1) % nlon).ravel()], axis=1)
    south = np.stack([flat[:-1].ravel(), flat[1:].ravel()], axis=1)
    return east[:, 0].astype(np.int64), east[:, 1].astype(np.int64), \
        south[:, 0].astype(np.int64), south[:, 1].astype(np.int64)


def components(n: int, u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Connected-component label per node, over the edges given."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    if u.size == 0:
        return np.arange(n, dtype=np.int64)
    data = np.ones(u.size, dtype=np.int8)
    g = coo_matrix((data, (u, v)), shape=(n, n)).tocsr()
    _, labels = connected_components(g, directed=False)
    return labels


def component_areas(labels: np.ndarray, area: np.ndarray,
                    qualifying: np.ndarray) -> np.ndarray:
    """Area of each component, counting only the qualifying nodes."""
    if not qualifying.any():
        return np.zeros(0)
    lab = labels[qualifying]
    _, inv = np.unique(lab, return_inverse=True)
    return np.bincount(inv, weights=area[qualifying])


def mesh_ocean_edges(export: Export, ocean: np.ndarray):
    """Undirected mesh edges with both ends in the ocean, and their sill depth.

    The sill of an edge is the SHALLOWER of its two regions, because that is what
    a flow between them has to cross. Nothing here is a grid quantity: this is
    the process-then-aggregate order's own graph.
    """
    off, lst = export.adjacency
    degree = np.diff(off.astype(np.int64))
    # int32 throughout: the mesh has under 2^31 regions and the half-edge list
    # runs to tens of millions, so the width is the difference between a
    # comfortable working set and a swapping one.
    src = np.repeat(np.arange(export.n_regions, dtype=np.int32), degree)
    dst = lst.astype(np.int32, copy=False)
    keep = src < dst
    np.logical_and(keep, ocean[src], out=keep)
    np.logical_and(keep, ocean[dst], out=keep)
    return src[keep], dst[keep]


def volume_arm(cell, ncell, area, depth, ocean, planet_area) -> dict:
    """CONTROL. Ocean volume is affine in depth, so the reduction is exact."""
    mean, ocean_area, covered = gridding.cell_mean(cell, ncell, area, depth, ocean)
    fine = float(np.sum(area[ocean] * depth[ocean]))
    coarse = float(np.sum(mean * ocean_area))
    scale = float(np.sum(area[ocean] * np.abs(depth[ocean])))
    return {
        "fine_volume_km3": fine,
        "coarse_volume_km3": coarse,
        "residual_relative": abs(coarse - fine) / scale if scale > 0 else 0.0,
        "tolerance": VOLUME_TOLERANCE,
        "passes": abs(coarse - fine) / scale <= VOLUME_TOLERANCE if scale > 0 else True,
        "covered_cells": int(covered.sum()),
        "ocean_area_fraction_of_planet": float(ocean_area.sum() / planet_area),
    }


def horizon_arm(cell, ncell, area, depth, ocean, planet_area,
                horizons, area_bar_km2) -> list[dict]:
    """The area on either side of a depth horizon, in both orders.

    `depth_only` isolates the threshold's own Jensen term: the cell keeps its
    true ocean area and only the test moves onto the cell mean. `binarised` adds
    the mask the atmosphere actually applies, so the pair separates the two
    mechanisms instead of reporting their sum.
    """
    mean, ocean_area, covered = gridding.cell_mean(cell, ncell, area, depth, ocean)
    total_area = gridding.cell_sum(cell, ncell, area)
    wet = np.zeros(ncell)
    np.divide(ocean_area, total_area, out=wet, where=total_area > 0)
    binary_ocean = wet >= 0.5

    # P3's bound: the ocean area a horizon can possibly move is the ocean area in
    # the cells the horizon CUTS. Both the shallowest and the deepest region of a
    # cell come from the DISTRIBUTION operator rather than a second reduction.
    quant, _ = gridding.cell_quantiles(cell, ncell, area, depth,
                                       probs=(0.0, 1.0), population=ocean)
    shallowest, deepest = quant[:, 0], quant[:, 1]

    out = []
    for h in horizons:
        fine = float(np.sum(area[ocean & (depth >= h)]))
        depth_only = float(np.sum(ocean_area[covered & (mean >= h)]))
        binarised = float(np.sum(total_area[binary_ocean & (mean >= h)]))
        straddling = covered & (shallowest < h) & (deepest >= h)
        bound = float(np.sum(ocean_area[straddling]))
        out.append({
            "horizon_km": h,
            "fine_area_km2": fine,
            "depth_only_area_km2": depth_only,
            "binarised_area_km2": binarised,
            "depth_only_gap_km2": depth_only - fine,
            "binarised_gap_km2": binarised - fine,
            "depth_only_gap_fraction_of_planet": (depth_only - fine) / planet_area,
            "binarised_gap_fraction_of_planet": (binarised - fine) / planet_area,
            "straddling_ocean_area_km2": bound,
            "gap_within_straddling_bound": abs(depth_only - fine) <= bound * (1 + 1e-9),
            "depth_only_gap_above_bar": abs(depth_only - fine) > area_bar_km2,
            "binarised_gap_above_bar": abs(binarised - fine) > area_bar_km2,
        })
    return out


def grid_connectivity(shape, mean, cell_area_grid, is_ocean, horizons,
                      planet_area, area_bar_km2) -> list[dict]:
    """Largest connected deep body per horizon on a grid, by area.

    Reported as the area NOT attached to the largest connected body, which is
    what a circulation cares about: a pool the terrain isolates and the grid
    attaches is a pool the model ventilates when the planet does not.
    """
    nlat, nlon = shape
    eu1, ev1, eu2, ev2 = grid_edges(nlat, nlon)
    u = np.concatenate([eu1, eu2])
    v = np.concatenate([ev1, ev2])
    out = []
    for h in horizons:
        deep = is_ocean & (mean >= h)
        keep = deep[u] & deep[v]
        labels = components(nlat * nlon, u[keep], v[keep])
        areas = component_areas(labels, cell_area_grid, deep)
        total = float(areas.sum())
        largest = float(areas.max()) if areas.size else 0.0
        out.append({
            "horizon_km": h,
            "deep_area_km2": total,
            "largest_body_km2": largest,
            "detached_area_km2": total - largest,
            "detached_fraction_of_planet": (total - largest) / planet_area,
            "bodies_above_area_bar": int((areas > area_bar_km2).sum()),
        })
    return out


def mesh_connectivity(n, area, depth, ocean, u, v, horizons,
                      planet_area, area_bar_km2) -> list[dict]:
    """The same spectrum on the mesh graph, which is the process-then-aggregate order."""
    sill = np.minimum(depth[u], depth[v])
    out = []
    for h in horizons:
        deep = ocean & (depth >= h)
        keep = sill >= h
        labels = components(n, u[keep], v[keep])
        areas = component_areas(labels, area, deep)
        total = float(areas.sum())
        largest = float(areas.max()) if areas.size else 0.0
        out.append({
            "horizon_km": h,
            "deep_area_km2": total,
            "largest_body_km2": largest,
            "detached_area_km2": total - largest,
            "detached_fraction_of_planet": (total - largest) / planet_area,
            "bodies_above_area_bar": int((areas > area_bar_km2).sum()),
        })
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", default=None,
                    help="build name; the configured one by default")
    ap.add_argument("--rungs", default="T21,T42,T85,T127,T170")
    ap.add_argument("--goldstein", default="36x36",
                    help="imaxXjmax for the OCN-3 candidate ocean grid")
    ap.add_argument("--igrid", type=int, default=gridding.GOLDSTEIN_EQUAL_AREA)
    ap.add_argument("--no-connectivity", action="store_true",
                    help="skip the mesh-graph arm, which is the expensive one")
    ap.add_argument("--out", type=Path, default=OUT_JSON)
    args = ap.parse_args()

    config = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    loop = yaml.safe_load(TRANSPORT_LOOP.read_text(encoding="utf-8"))
    area_criterion = _criterion(loop, "sea_ice_area")
    heat_criterion = _criterion(loop, "heat_transport_change")

    # The mesh lives in the T42 export's `raw/`, and the build is the configured
    # one unless the caller names another. `CLAUDE.md` rule 5: the pointer is
    # deliberate, and `builds.mesh_export` is the door that resolves it.
    root = (builds.mesh_export(config) if args.build is None
            else builds.mesh_export_at(ROOT / "source" / args.build))
    export = Export(root)

    radius_km = export.radius_km
    planet_area_km2 = 4.0 * np.pi * radius_km ** 2
    area_bar_km2 = float(area_criterion["threshold"]) * planet_area_km2

    ocean = export.surface_class == OCEAN
    depth = np.maximum(-export.elevation_km.astype(np.float64), 0.0)
    area = export.cell_area.astype(np.float64)
    n = export.n_regions
    ocean_area_km2 = float(area[ocean].sum())

    supports: list[tuple[str, np.ndarray, int, tuple[int, int], np.ndarray]] = []
    checks: list[dict] = []
    for rung in [r.strip().upper() for r in args.rungs.split(",") if r.strip()]:
        nlat, nlon, _ = rungs.geometry(rung)
        grid_dir = root.parent / f"exoplasim-{rung}"
        cell, gnlat, gnlon = gridding.region_cells(export, grid_dir)
        if (gnlat, gnlon) != (nlat, nlon):
            raise SystemExit(f"{rung}: the export grid is {gnlat}x{gnlon}, "
                             f"the ladder says {nlat}x{nlon}")
        spec = gridding.gaussian_grid(nlat, nlon, name=f"exoplasim-{rung}")
        supports.append((f"exoplasim-{rung}", cell, nlat * nlon, (nlat, nlon),
                         spec.cell_area(radius_km * 1000.0).ravel() / 1.0e6))

    ni, nj = (int(v) for v in args.goldstein.lower().split("x"))
    ocn_spec = gridding.goldstein_grid(ni, nj, igrid=args.igrid)
    gcell = spec_cells(ocn_spec, export.lat, export.lon)
    supports.append((ocn_spec.name, gcell, ocn_spec.ncell, (nj, ni),
                     ocn_spec.cell_area(radius_km * 1000.0).ravel() / 1.0e6))

    # The binner is the one coordinate source, and this is the identity that says
    # so: on a spec whose columns ARE the export's, it must reproduce
    # `gridding.column` for every region, exactly. A wrap that drifted by a cell
    # would show here rather than in a result nobody can check.
    probe = gridding.gaussian_grid(64, 128, name="column-probe")
    mine = spec_cells(probe, export.lat, export.lon) % 128
    theirs = gridding.column(export.lon, 128)
    checks.append({
        "check": "the spec binner reproduces gridding.column on the export's own columns",
        "mismatched_regions": int(np.sum(mine != theirs)),
        "tolerance": 0,
        "passes": bool(np.all(mine == theirs)),
    })

    mesh_u = mesh_v = None
    if not args.no_connectivity:
        mesh_u, mesh_v = mesh_ocean_edges(export, ocean)

    report = {
        "build": root.parent.name,
        "terrain_hash": export.provenance().get("terrain_hash"),
        "regions": int(n),
        "radius_km": radius_km,
        "planet_area_km2": planet_area_km2,
        "ocean_area_km2": ocean_area_km2,
        "ocean_area_fraction": ocean_area_km2 / planet_area_km2,
        "horizons_km": list(HORIZONS_KM),
        "bars": {
            "volume_relative": {
                "value": VOLUME_TOLERANCE,
                "from": "lib/remap.py:CLOSURE_TOLERANCE, the float64 round-off floor",
            },
            "area_fraction_of_planet": {
                "value": float(area_criterion["threshold"]),
                "km2": area_bar_km2,
                "from": "ocean/config/transport_loop.yaml, criterion sea_ice_area",
            },
            "heat_transport_w_m2": {
                "value": float(heat_criterion["threshold"]),
                "from": "ocean/config/transport_loop.yaml, criterion heat_transport_change",
            },
        },
        "checks": checks,
        "supports": {},
    }

    for label, cell, ncell, shape, cell_area_grid in supports:
        mean, ocean_area, covered = gridding.cell_mean(cell, ncell, area, depth, ocean)
        total_area = gridding.cell_sum(cell, ncell, area)
        wet = np.zeros(ncell)
        np.divide(ocean_area, total_area, out=wet, where=total_area > 0)
        binary_ocean = wet >= 0.5

        entry = {
            "shape": list(shape),
            "cells": int(ncell),
            "cells_with_ocean": int(covered.sum()),
            "cells_binarised_ocean": int(binary_ocean.sum()),
            "volume": volume_arm(cell, ncell, area, depth, ocean, planet_area_km2),
            "horizon_area": horizon_arm(cell, ncell, area, depth, ocean,
                                        planet_area_km2, HORIZONS_KM, area_bar_km2),
        }

        # P2. A uniform ocean reduces to itself, so the horizon gap is exactly
        # zero at every horizon and every support. This is `cell_mean`'s own
        # constant-preservation invariant reached through this arm.
        flat = np.where(ocean, 1.0, 0.0)
        flat_mean, _, flat_cov = gridding.cell_mean(cell, ncell, area, flat, ocean)
        worst = 0.0
        for h in (0.5, 1.0, 1.5):
            fine = float(np.sum(area[ocean & (flat >= h)]))
            coarse = float(np.sum(ocean_area[flat_cov & (flat_mean >= h)]))
            worst = max(worst, abs(coarse - fine))
        entry["uniform_control"] = {
            "worst_absolute_gap_km2": worst,
            "tolerance_km2": 0.0,
            "passes": worst == 0.0,
        }

        if not args.no_connectivity:
            entry["connectivity_coarse"] = grid_connectivity(
                shape, mean, total_area, binary_ocean, HORIZONS_KM,
                planet_area_km2, area_bar_km2)
        report["supports"][label] = entry

    if not args.no_connectivity:
        report["connectivity_fine"] = mesh_connectivity(
            n, area, depth, ocean, mesh_u, mesh_v, HORIZONS_KM,
            planet_area_km2, area_bar_km2)

    # The capacity bar in its own units. A passage carrying an ocean heat
    # transport H whose cross-section the reduction moves by a fraction f moves
    # the heat by f*H, and the loop's controlling scalar is reached at
    # f*H = threshold * ocean area. H is not knowable without an ocean run, so
    # what is published is the critical fraction against it.
    ocean_area_m2 = ocean_area_km2 * 1.0e6
    material_w = float(heat_criterion["threshold"]) * ocean_area_m2
    report["capacity_bar"] = {
        "ocean_area_m2": ocean_area_m2,
        "material_heat_transport_w": material_w,
        "material_heat_transport_pw": material_w / 1.0e15,
        "critical_fraction_by_passage_transport_pw": {
            f"{h:g}": material_w / (h * 1.0e15) for h in (0.1, 0.25, 0.5, 1.0, 2.0)
        },
        "decision_procedure": (
            "Take the passage's own heat transport from the first accepted "
            "ocean run and compare the measured fractional capacity error "
            "against material_heat_transport_w divided by it. Above one the "
            "passage needs its geometry carried at a finer support; at or "
            "below one the cell-mean bathymetry is admissible there."),
        "sea_water_density_kg_m3": SEA_WATER_DENSITY_KG_M3,
    }
    report["provenance"] = provenance.config_stamp(
        config, "analysis/ocean_support_reduction.py")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    failed = [c for c in checks if not c["passes"]]
    for label, entry in report["supports"].items():
        v = entry["volume"]
        mark = "  ok  " if v["passes"] else " FAIL "
        print(f"[{mark}] {label}: volume residual {v['residual_relative']:.3g} "
              f"<= {VOLUME_TOLERANCE:.0e}")
        u = entry["uniform_control"]
        mark = "  ok  " if u["passes"] else " FAIL "
        print(f"[{mark}] {label}: uniform-ocean horizon gap "
              f"{u['worst_absolute_gap_km2']:.3g} km2, must be exactly zero")
        if not v["passes"]:
            failed.append(v)
        if not u["passes"]:
            failed.append(u)
    for c in checks:
        mark = "  ok  " if c["passes"] else " FAIL "
        print(f"[{mark}] {c['check']}")
    print(f"\n{paths_lib.rel(args.out)}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
