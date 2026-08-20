#!/usr/bin/env python3
"""Build the climate-independent hydrography products.

Everything here depends only on geometry, so it can be built before any climate
exists. What it deliberately does not do is decide lake levels: that is a water
balance, it needs precipitation and evaporation, and it belongs to
`lake_balance.py` once ExoPlaSim has run.

Products, all under `hydrography/data/<build>/`:

  regions.nc      per mesh region: terminal, receiver and filled surface
  basins.nc       per preserved basin: final-terrain hypsometry, spill, catchment
  coupling_<grid>.nc  sparse basin-by-grid-cell catchment areas, the interface a
                  climate model integrates precipitation over
  hydrography_report.json  diagnostics and provenance
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from netCDF4 import Dataset
import numpy as np

# _paths first: it is what puts the project's lib/ on the path, so anything
# imported from there has to come after it.
from _paths import ANALYSIS, DATA
from builds import build_root, grid_export
import drainage as dr
import gridding
from orogen import Export, LAND, OCEAN

N_LEVELS = 128  # samples per basin hypsometric curve


def hypsometry(export: Export, drn: dr.Drainage, spill: np.ndarray, n_basins: int):
    """Final-terrain level/area/volume curves, one per basin.

    Derived from the flood rather than from the catalogue. The catalogue's
    `hypsometry` is measured on the natural, pre-conditioning terrain and
    overstates capacity badly once erosion has worked on the rim: for the first
    basin it gives 53,968 km2 flooded at spill against 10,853 km2 on the
    finished surface.

    The flood already visited every region in ascending order of the level a
    droplet must reach to get there, which is exactly the hypsometric ordering,
    so the curve is a cumulative sum over the basin's members below its spill.
    """
    area = export.cell_area.astype(np.float64)
    term, filled = drn.terminal, drn.filled_km

    levels = np.zeros((n_basins, N_LEVELS))
    areas = np.zeros((n_basins, N_LEVELS))
    volumes = np.zeros((n_basins, N_LEVELS))

    # Group by basin first, then order by level within each group. Sorting by
    # level alone and grouping with searchsorted silently returns nonsense,
    # because searchsorted needs the key it bisects to be the sorted one.
    sel = np.flatnonzero(term >= 0)
    order = sel[np.lexsort((filled[sel], term[sel]))]
    bounds = np.searchsorted(term[order], np.arange(n_basins + 1))

    for b in range(n_basins):
        members = order[bounds[b]:bounds[b + 1]]
        if members.size == 0:
            continue
        lv, ar = filled[members], area[members]
        keep = lv <= spill[b]
        lv, ar = lv[keep], ar[keep]
        if lv.size == 0:                      # sink already at or above spill
            levels[b] = spill[b]
            continue
        # Evaluate both curves exactly at the sample levels. Flooded area is a
        # step function of level, so interpolating it between samples is wrong
        # by up to 80% in basins that spread suddenly across a flat, and those
        # are precisely the basins where lake extent is most sensitive.
        #   A(L) = sum of a_i where f_i <= L
        #   V(L) = sum of (L - f_i) a_i = L*A(L) - sum(f_i a_i)
        cum_area = np.cumsum(ar)
        cum_fa = np.cumsum(lv * ar)
        grid = np.linspace(lv[0], spill[b], N_LEVELS)
        idx = np.searchsorted(lv, grid, side="right")
        A = np.where(idx > 0, cum_area[np.maximum(idx - 1, 0)], 0.0)
        FA = np.where(idx > 0, cum_fa[np.maximum(idx - 1, 0)], 0.0)
        levels[b] = grid
        areas[b] = A
        volumes[b] = grid * A - FA
    return levels, areas, volumes


def couple_to_grid(export: Export, drn: dr.Drainage, grid_dir, n_basins: int):
    """Sparse basin-by-grid-cell catchment areas.

    This is the interface to a climate model: to get net inflow to a basin you
    integrate (P - E) over its catchment, and the catchment crosses grid cells.
    Regions are assigned to the cell containing their centre, matching the
    exporter's own categorical resampling rule.
    """
    # The cell index comes from `lib/gridding.py`, which owns the convention.
    # This function used to derive its own copy of the same two expressions,
    # which is how the export and the model came to be reconciled by longitude
    # LABEL twice: with the arithmetic written out in four places there was no
    # single thing to be right, only four things to agree.
    lat, lon, gm = gridding.grid_geometry(grid_dir)
    nlat, nlon = lat.size, lon.size
    cell, _, _ = gridding.region_cells(export, grid_dir)

    land = export.surface_class == LAND
    sel = land & (drn.terminal >= 0)
    key = drn.terminal[sel].astype(np.int64) * (nlat * nlon) + cell[sel]
    uniq, inv = np.unique(key, return_inverse=True)
    acc = np.zeros(uniq.size)
    np.add.at(acc, inv, export.cell_area[sel].astype(np.float64))
    return {
        "cell_lat": lat, "cell_lon": lon,
        "n_lat": nlat, "n_lon": nlon,
        "basin": (uniq // (nlat * nlon)).astype(np.int32),
        "cell": (uniq % (nlat * nlon)).astype(np.int32),
        "area_km2": acc,
        "grid_name": grid_dir.name,
        "truncation": gm.get("truncation"),
    }


def cross_check_against_export(export: Export, drn: dr.Drainage) -> dict:
    """Compare our routing against the export's own, and record the agreement.

    The export now routes drainage itself, correctly, and we could in principle
    read `drainage_terminal` instead of flooding. We do not, for two reasons:
    the flood also yields the filled surface that the hypsometry curves are built
    from, and an independent implementation is what caught two routing bugs in
    the export. Keeping both and recording their agreement turns that into a
    standing regression check rather than a one-off comparison.
    """
    sc, area, dt = export.surface_class, export.cell_area.astype(np.float64), export.drainage_terminal
    land = sc == LAND
    sinks = np.array([b.sink for b in export.basins])
    lut = np.full(export.n_regions, -1, np.int32)
    lut[sinks] = np.arange(len(sinks))
    theirs = np.full(export.n_regions, dr.TERMINAL_OCEAN, np.int32)
    is_sink = np.isin(dt, sinks)
    theirs[is_sink] = lut[dt[is_sink]]

    mine = drn.catchment_areas(export, len(export.basins))
    declared = np.array([b.catchment_area_km2 for b in export.basins])
    ratio = mine / np.maximum(declared, 1e-9)

    # The export now publishes its own endorheic fraction with the denominator
    # named. Assert against the published value rather than recomputing it, so a
    # future change of denominator on either side shows up here instead of being
    # silently absorbed.
    dc = export.manifest["basins"].get("drainageConsistency", {})
    published = dc.get("fractionOfLand")
    ours_by_export_routing = float(
        area[land & (theirs >= 0)].sum() / area[land].sum())
    if published is not None and abs(published - ours_by_export_routing) > 1e-6:
        raise RuntimeError(
            f"export publishes fractionOfLand={published} but its own "
            f"drainage_terminal gives {ours_by_export_routing}; the denominator "
            f"is {dc.get('landDenominator')!r}"
        )

    return {
        "region_agreement": float((theirs[land] == drn.terminal[land]).mean()),
        "regions_differing": int((theirs[land] != drn.terminal[land]).sum()),
        "catchment_ratio_median": float(np.median(ratio)),
        "catchment_ratio_p5": float(np.percentile(ratio, 5)),
        "catchment_ratio_p95": float(np.percentile(ratio, 95)),
        "endorheic_land_fraction_ours": float(
            area[land & (drn.terminal >= 0)].sum() / area[land].sum()),
        "endorheic_land_fraction_export": ours_by_export_routing,
        "endorheic_land_fraction_export_published": published,
        "export_land_denominator": dc.get("landDenominator"),
        "note": ("Land fractions are normalised by surface_class land area. The "
                 "export now publishes its own fraction with the denominator "
                 "named, and this build asserts the two agree."),
    }


def river_mouths(export: Export, drn: dr.Drainage, top: int = 200):
    """Largest discharges reaching the world ocean, on the resolved network."""
    off, adj = export.adjacency
    sc, acc = export.surface_class, export.flow_accumulation
    n = export.n_regions
    src = np.repeat(np.arange(n, dtype=np.int64), np.diff(off).astype(np.int64))
    coastal = np.zeros(n, dtype=bool)
    np.logical_or.at(coastal, src, sc[adj] == OCEAN)
    mouth = coastal & (sc == LAND) & (drn.terminal == dr.TERMINAL_OCEAN)
    idx = np.flatnonzero(mouth)
    idx = idx[np.argsort(acc[idx])[::-1][:top]]
    return [{"region": int(i), "lat": float(export.lat[i]), "lon": float(export.lon[i]),
             "discharge_area_km2": float(acc[i])} for i in idx]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--export", type=str, default=None, help="Orogen export with raw/")
    ap.add_argument("--grids", nargs="*", default=["exoplasim-T42", "exoplasim-T85"])
    ap.add_argument("--output", type=str, default=None)
    args = ap.parse_args()

    ANALYSIS.mkdir(parents=True, exist_ok=True)

    ex = Export(args.export)
    n_basins = len(ex.basins)
    terrain_build = ex.provenance()["terrain_build"]
    print(f"export: {terrain_build}")

    # Output is namespaced by the build being PROCESSED, not by config's active
    # build and not by the flat data/. This script is the one that creates a
    # per-build directory, so it takes the name from the export in hand -- if it
    # took it from config, building a non-active build would overwrite the
    # active one's products under the active one's name.
    if args.output is not None:
        out = Path(args.output)
    else:
        if terrain_build == "unrecognised":
            raise SystemExit(
                "export terrain is not a registered build, so its products have "
                "no name to be filed under; register it in lib/orogen.py or pass "
                "--output explicitly")
        out = DATA / terrain_build
    out.mkdir(parents=True, exist_ok=True)
    print(f"{ex.n_regions:,} regions, {n_basins} preserved basins")

    print("resolving drainage (priority flood)...")
    drn = dr.resolve(ex)
    land = ex.surface_class == LAND
    if (land & (drn.terminal == dr.TERMINAL_NONE)).any():
        raise RuntimeError("priority flood left land regions unresolved")

    agreement = cross_check_against_export(ex, drn)
    print(f"agreement with the export's own drainage_terminal: "
          f"{agreement['region_agreement']*100:.3f}% of land, "
          f"catchment ratio median {agreement['catchment_ratio_median']:.4f}")

    spill, spill_target, spill_region, spill_exit, merged = dr.spill_levels(
        ex, drn, n_basins)
    if not np.isfinite(spill).all():
        raise RuntimeError(
            f"{int((~np.isfinite(spill)).sum())} basins have no finite spill level; "
            "the hypsometry would be meaningless"
        )
    catch = drn.catchment_areas(ex, n_basins)
    print("building final-terrain hypsometry...")
    levels, areas, volumes = hypsometry(ex, drn, spill, n_basins)

    area = ex.cell_area.astype(np.float64)
    total = area.sum()
    endorheic_land = area[land & (drn.terminal >= 0)].sum()

    # -- regions.nc ---------------------------------------------------------
    with Dataset(out / "regions.nc", "w") as ds:
        ds.createDimension("region", ex.n_regions)
        ds.title = "Resolved drainage on the World Orogen mesh"
        ds.terrain_hash = ex.terrain_hash
        v = ds.createVariable("terminal", "i4", ("region",), zlib=True)
        v.long_name = "preserved basin index, -1 world ocean, -2 not land"
        v[:] = drn.terminal
        v = ds.createVariable("receiver", "i4", ("region",), zlib=True)
        v.long_name = ("region this one drains into, -1 at the world ocean or a "
                       "basin sink; the priority flood's discovery pointer, which "
                       "routes across filled flats where steepest descent cannot")
        v[:] = drn.receiver
        v = ds.createVariable("filled_km", "f4", ("region",), zlib=True)
        v.long_name = "depression-filled surface elevation"
        v.units = "km"
        v[:] = drn.filled_km

    # -- basins.nc ----------------------------------------------------------
    with Dataset(out / "basins.nc", "w") as ds:
        ds.createDimension("basin", n_basins)
        ds.createDimension("level", N_LEVELS)
        ds.title = "Preserved closed basins, hypsometry on the finished terrain"
        ds.terrain_hash = ex.terrain_hash
        ds.note = ("Curves are recomputed from the resolved drainage, not taken "
                   "from manifest.basins.preserved[].hypsometry, which is "
                   "measured on the natural pre-conditioning terrain.")
        for name, data, units, desc in [
            ("sink", np.array([b.sink for b in ex.basins], "i4"), "1", "sink region index"),
            ("sink_elevation_km", np.array([b.sink_elevation_km for b in ex.basins]), "km", "sink elevation"),
            ("spill_km", spill, "km", "spill level on the finished terrain"),
            ("catchment_km2", catch, "km2", "land area draining to this basin"),
            ("capacity_km3", volumes[:, -1], "km3", "volume held at spill"),
            ("area_at_spill_km2", areas[:, -1], "km2", "flooded area at spill"),
            ("spill_target", spill_target, "1", "basin index this overflows into, -1 world ocean"),
            ("spill_region", spill_region, "1",
             "region on the basin side of the saddle it overflows at, -1 if it has no exit"),
            ("spill_exit_region", spill_exit, "1",
             "region on the far side of that saddle, where the outflow river starts; "
             "may be ocean, in which case the overflow enters the sea directly"),
            ("critical_aridity_index", catch / np.maximum(areas[:, -1], 1e-9) - 1.0, "1",
             "basin overflows, and so should carve its outlet, wherever (E-P)/runoff "
             "over its catchment falls below this"),
        ]:
            v = ds.createVariable(name, data.dtype if data.dtype.kind == "i" else "f8",
                                  ("basin",), zlib=True)
            v.units, v.long_name = units, desc
            v[:] = data
        for name, data, units in [("level_km", levels, "km"),
                                  ("flooded_area_km2", areas, "km2"),
                                  ("volume_km3", volumes, "km3")]:
            v = ds.createVariable(name, "f8", ("basin", "level"), zlib=True)
            v.units = units
            v[:] = data
        ids = ds.createVariable("basin_id", str, ("basin",))
        for i, b in enumerate(ex.basins):
            ids[i] = b.id

    # -- coupling matrices --------------------------------------------------
    couplings = {}
    for g in args.grids:
        gd = build_root() / g
        if not (gd / "manifest.json").is_file():
            print(f"  skipping {g}: not present")
            continue
        c = couple_to_grid(ex, drn, gd, n_basins)
        couplings[g] = {"pairs": int(c["basin"].size), "n_lat": c["n_lat"], "n_lon": c["n_lon"]}
        with Dataset(out / f"coupling_{g}.nc", "w") as ds:
            ds.createDimension("pair", c["basin"].size)
            ds.title = f"Basin catchment area per {g} grid cell (sparse COO)"
            ds.terrain_hash = ex.terrain_hash
            ds.n_lat, ds.n_lon = c["n_lat"], c["n_lon"]
            ds.note = ("cell index is row * n_lon + col, rows north to south, "
                       "on the columns lib/gridding.py:column() assigns. Read "
                       "climate fields by index, never by longitude label; "
                       "gridding.require_index_alignment is the check.")
            ds.createDimension("cell_lat", c["n_lat"])
            ds.createDimension("cell_lon", c["n_lon"])
            for nm, dat in [("cell_lat", c["cell_lat"]), ("cell_lon", c["cell_lon"])]:
                v = ds.createVariable(nm, "f8", (nm,))
                v.units = "degrees_north" if nm.endswith("lat") else "degrees_east"
                v.long_name = "coordinate the cell index was built on"
                v[:] = dat
            for nm, dat, ty, un in [("basin", c["basin"], "i4", "1"),
                                    ("cell", c["cell"], "i4", "1"),
                                    ("area_km2", c["area_km2"], "f8", "km2")]:
                v = ds.createVariable(nm, ty, ("pair",), zlib=True)
                v.units = un
                v[:] = dat
        print(f"  coupled to {g}: {c['basin'].size:,} basin-cell pairs")

    # -- report -------------------------------------------------------------
    report = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "source": ex.provenance(),
        "drainage": {
            "land_fraction_of_planet": float(area[land].sum() / total),
            "endorheic_fraction_of_land": float(endorheic_land / area[land].sum()),
            "ocean_draining_fraction_of_land": float(1 - endorheic_land / area[land].sum()),
            "unresolved_regions": 0,
            "pits_filled": int((drn.filled_km > ex.elevation_km)[land].sum()),
            "median_fill_depth_m": float(np.median(
                ((drn.filled_km - ex.elevation_km)[land])[
                    (drn.filled_km - ex.elevation_km)[land] > 0]) * 1000)
            if (drn.filled_km > ex.elevation_km)[land].any() else 0.0,
        },
        "basins": {
            "count": n_basins,
            "total_capacity_km3": float(volumes[:, -1].sum()),
            "total_area_at_spill_km2": float(areas[:, -1].sum()),
            "area_at_spill_fraction_of_planet": float(areas[:, -1].sum() / total),
            "total_catchment_km2": float(catch.sum()),
            "catchment_to_spill_area_ratio_median": float(
                np.median(catch / np.maximum(areas[:, -1], 1e-9))),
            "merged_basins": int(merged),
            "spills_to_ocean": int((spill_target < 0).sum()),
            "spills_into_another_basin": int((spill_target >= 0).sum()),
            "critical_aridity_index_percentiles": {
                str(q): float(v) for q, v in zip(
                    (5, 25, 50, 75, 95),
                    np.percentile(catch / np.maximum(areas[:, -1], 1e-9) - 1.0,
                                  [5, 25, 50, 75, 95]))},
            # Both sides are on the relief-scaled vertical from the build that
            # fixed the catalogue's ...Km fields onward. Before that the
            # denominator was unscaled and this ratio drifted with gravity: 4%
            # off at 10.1989 m/s2 and 31% at 12.81, in a manifest that stated
            # one unit for both. Builds predating the fix cannot be corrected by
            # multiplying through, because the scaling follows elevation_km's
            # own rule and does not apply below sea level.
            "capacity_vs_natural_catalogue": float(
                volumes[:, -1].sum()
                / sum(b.natural_volume_km3 for b in ex.basins)),
        },
        "coupling": couplings,
        "cross_check_vs_export_routing": agreement,
        "river_mouths_top": river_mouths(ex, drn, top=50),
        # Sill depth decides whether a marginal sea exchanges freely or turns
        # evaporitic or stratified. The full list stays in the source manifest;
        # only the largest are echoed here to keep the report readable.
        "marginal_seas": {
            "count": len(ex.manifest["hydrology"].get("oceanBasins") or []),
            "largest": sorted(
                (ex.manifest["hydrology"].get("oceanBasins") or []),
                key=lambda s: -(s.get("areaKm2") or 0.0))[:25],
            "source": "manifest.hydrology.oceanBasins",
        },
    }
    (out / "hydrography_report.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")

    d = report["drainage"]; b = report["basins"]
    print(f"\nland {d['land_fraction_of_planet']*100:.2f}% of planet")
    print(f"  endorheic share of land : {d['endorheic_fraction_of_land']*100:.1f}%")
    print(f"  pits filled             : {d['pits_filled']:,} regions, "
          f"median {d['median_fill_depth_m']:.1f} m")
    print(f"basins: capacity {b['total_capacity_km3']:,.0f} km3, "
          f"area at spill {b['area_at_spill_fraction_of_planet']*100:.2f}% of planet")
    print(f"  that capacity is {b['capacity_vs_natural_catalogue']*100:.1f}% of the "
          f"catalogue's natural-terrain figure")
    print(f"\nwrote products to {out}")


if __name__ == "__main__":
    main()
