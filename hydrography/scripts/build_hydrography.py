#!/usr/bin/env python3
"""Build the climate-independent hydrography products.

Everything here depends only on geometry, so it can be built before any climate
exists. What it deliberately does not do is decide lake levels: that is a water
balance, it needs precipitation and evaporation, and it belongs to
`lake_balance.py` once ExoPlaSim has run.

Products, all under `hydrography/data/`:

  regions.nc      per mesh region: resolved drainage terminal and filled surface
  basins.nc       per preserved basin: final-terrain hypsometry, spill, catchment
  coupling_<grid>.nc  sparse basin-by-grid-cell catchment areas, the interface a
                  climate model integrates precipitation over
  hydrography_report.json  diagnostics and provenance
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json

from netCDF4 import Dataset
import numpy as np

from _paths import ANALYSIS, DATA, SOURCE
import drainage as dr
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
    gm = json.loads((grid_dir / "manifest.json").read_text(encoding="utf-8"))["grid"]
    lat = np.fromfile(grid_dir / gm["coords"]["lat"]["path"], dtype="float64")
    lon = np.fromfile(grid_dir / gm["coords"]["lon"]["path"], dtype="float64")
    nlat, nlon = lat.size, lon.size

    # Latitude rows run north to south and may be Gaussian, so bin on the
    # midpoints between row centres rather than assuming uniform spacing.
    edges = np.empty(nlat + 1)
    edges[1:-1] = 0.5 * (lat[:-1] + lat[1:])
    edges[0], edges[-1] = 90.0, -90.0
    rlat, rlon = export.lat, export.lon
    row = np.clip(np.searchsorted(-edges, -rlat, side="right") - 1, 0, nlat - 1)
    col = np.clip(((rlon + 180.0) / 360.0 * nlon).astype(np.int64), 0, nlon - 1)
    cell = row * nlon + col

    land = export.surface_class == LAND
    sel = land & (drn.terminal >= 0)
    key = drn.terminal[sel].astype(np.int64) * (nlat * nlon) + cell[sel]
    uniq, inv = np.unique(key, return_inverse=True)
    acc = np.zeros(uniq.size)
    np.add.at(acc, inv, export.cell_area[sel].astype(np.float64))
    return {
        "n_lat": nlat, "n_lon": nlon,
        "basin": (uniq // (nlat * nlon)).astype(np.int32),
        "cell": (uniq % (nlat * nlon)).astype(np.int32),
        "area_km2": acc,
        "grid_name": grid_dir.name,
        "truncation": gm.get("truncation"),
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

    out = DATA if args.output is None else __import__("pathlib").Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    ANALYSIS.mkdir(parents=True, exist_ok=True)

    ex = Export(args.export)
    n_basins = len(ex.basins)
    print(f"export: {ex.provenance()['terrain_build']}")
    print(f"{ex.n_regions:,} regions, {n_basins} preserved basins")

    print("resolving drainage (priority flood)...")
    drn = dr.resolve(ex)
    land = ex.surface_class == LAND
    if (land & (drn.terminal == dr.TERMINAL_NONE)).any():
        raise RuntimeError("priority flood left land regions unresolved")

    spill, spill_target, merged = dr.spill_levels(ex, drn, n_basins)
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
        gd = SOURCE / g
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
            ds.note = "cell index is row * n_lon + col, rows north to south"
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
            "capacity_vs_natural_catalogue": float(
                volumes[:, -1].sum()
                / sum(b.natural_volume_km3 for b in ex.basins)),
        },
        "coupling": couplings,
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
