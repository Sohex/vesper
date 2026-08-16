"""Feed the pedology component's soil water capacity back to ExoPlaSim.

This closes the loop the other way. ExoPlaSim's land hydrology is a single
bucket per cell of depth `dwmax`, surface code 0229, and its runoff is literally
that bucket overflowing:

    landmod.f90:1098   drunoff(:) = AMAX1(0., dwatc(:) - dwmax(:)) / deltsec

`configure()` clears every surface field when handed a landmap, so `dwmax` falls
back to the uniform namelist default of `wsmax = 0.5 m` everywhere. That default
is why `model.uniform_land_surface` exists, and it is the most likely reason the
0.96 climatology reports a land runoff ratio of 2.8% where Earth's land manages
about 35%.

`pedology/` computes the real capacity from texture and regolith depth, which is
exactly the field being defaulted. Supplying it makes runoff a property of the
soil rather than of a namelist, and the soil is in turn a product of the climate
and the biosphere. That is the third loop in this pipeline.

    python exoplasim/scripts/build_surface_soil_water.py

**This changes climate results.** It is off unless `model.soil_water_source` is
set to `pedology` in `config/planet.yaml`, and that key is deliberately absent by
default: adding it moves `config_sha256` and blocks resumption of any run in
flight. Add it when re-baselining, not before.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

from _paths import CONFIG, INPUTS, PROJECT_ROOT
from convert_orogen import write_sra

SOIL_WATER_CODE = 229

# ExoPlaSim's own default, landmod.f90 `WSMAX_EARTH`. Non-land cells keep it,
# since dwmax is meaningless over ocean but must still be a sane number.
EXOPLASIM_DEFAULT_WSMAX_M = 0.5

# Coordinate rounding shared with pedology/scripts/build_soil.py.
COORD_DECIMALS = 4


def read_soil_map(path: Path) -> dict[tuple[float, float], float]:
    """Plant-available water capacity in mm, keyed by rounded coordinates."""
    lines = path.read_text().splitlines()
    header = lines[0].split()
    if "awc" not in header:
        raise SystemExit(
            f"{path} has no awc column. Rebuild it with "
            f"pedology/scripts/build_soil.py.")
    column = header.index("awc")
    capacity = {}
    for line in lines[1:]:
        parts = line.split()
        capacity[(round(float(parts[0]), COORD_DECIMALS),
                  round(float(parts[1]), COORD_DECIMALS))] = float(parts[column])
    return capacity


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--soil-map", type=Path,
                        default=None)
    parser.add_argument("--climatology", type=Path,
                        default=(PROJECT_ROOT / "exoplasim" / "analysis"
                                 / "climatology_s096"
                                 / "baseline_regular_climatology.nc"),
                        help="supplies the grid and the land mask, so they match "
                             "the soil map exactly")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--lakes", type=Path, default=None,
                        help="surface_water.nc; raises dwmax on the lake fraction "
                             "of each cell so those cells reach the wetness "
                             "ceiling and evaporate at the potential rate")
    parser.add_argument("--lake-dwmax-m", type=float, default=None,
                        help="bucket depth on the lake fraction; defaults to "
                             "model.lake_dwmax_m, else 0.2 m")
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    model = config["model"]
    if args.soil_map is None:
        import sys as _sys
        _sys.path.insert(0, str(PROJECT_ROOT / "lib"))
        import builds as _b
        args.soil_map = _b.soilmap(config)
    nlat, nlon = int(model["latitudes"]), int(model["longitudes"])
    resolution = str(model["resolution"]).upper()

    if not args.soil_map.is_file():
        raise SystemExit(
            f"{args.soil_map} does not exist. Run pedology/scripts/build_soil.py.")
    capacity_mm = read_soil_map(args.soil_map)

    # Grid and land mask come from the climatology, which is the boundary land
    # mask as the model itself saw it, and is the same grid pedology wrote its
    # soil map on. Deriving either independently is how the first version of
    # this script matched zero cells out of 4106.
    import netCDF4 as nc

    with nc.Dataset(args.climatology) as data:
        lat = np.asarray(data["lat"][:], dtype=float)
        lon = np.asarray(data["lon"][:], dtype=float)
        land = np.asarray(data["lsm"][0], dtype=float) > 0.5
    if (len(lat), len(lon)) != (nlat, nlon):
        raise SystemExit(
            f"climatology grid is {len(lat)}x{len(lon)} but config says "
            f"{nlat}x{nlon}")

    lon_signed = np.round(np.where(lon > 180.0, lon - 360.0, lon), COORD_DECIMALS)
    lat_rounded = np.round(lat, COORD_DECIMALS)

    field = np.full((nlat, nlon), EXOPLASIM_DEFAULT_WSMAX_M)
    matched = 0
    unmatched_land = 0
    for j in range(nlat):
        for i in range(nlon):
            if not land[j, i]:
                continue
            value = capacity_mm.get((float(lon_signed[i]), float(lat_rounded[j])))
            if value is None:
                unmatched_land += 1
                continue
            field[j, i] = value / 1000.0   # mm -> m, the units dwmax is in
            matched += 1

    # Lakes, as an area-weighted bucket depth. `drhs` reaches 1 once soil water
    # exceeds 40% of dwmax (landmod.f90:52-53), so a SHALLOWER bucket saturates
    # the wetness factor on less water and evaporates at the potential rate,
    # where a deeper one needs proportionally more water to get there. The note
    # in lake-representation.md reads the other way, "large and full"; large is
    # right for storage and wrong for wetness, and on this planet the lakes sit
    # in the arid cells where the water to fill a large bucket is not available.
    #
    # WHAT THIS CANNOT DO, verified in the source rather than assumed. A lake
    # here is fed by its catchment, and that water never reaches the evaporating
    # bucket. `dwatc` gains only from local precipitation minus evaporation
    # (landmod.f90:1097); routed river water accumulates into `driver`
    # (landmod.f90:1390), a separate store that discharges to the ocean. So a
    # cell's annual evaporation is capped by its own precipitation plus storage
    # however dwmax is set, and lake evaporation stays underestimated. What this
    # buys is the seasonal partition -- a lake cell that stays at potential rate
    # while it has water, rather than being moisture-limited from the first dry
    # day. Sustained open-water evaporation needs the mask flip in
    # notes/lake-representation.md, and only for the few resolvable lakes.
    lake_report = None
    if args.lakes is not None:
        import sys
        sys.path.insert(0, str(PROJECT_ROOT / "lib"))
        from gridding import land_fraction_of_class
        from orogen import Export
        from builds import grid_export, mesh_export
        mesh = Export(mesh_export(config))
        with nc.Dataset(args.lakes) as lds:
            lake = np.asarray(lds["lake"][:]).astype(bool)
            lake_terrain = getattr(lds, "terrain_hash", None)
        if lake_terrain and lake_terrain != mesh.terrain_hash:
            raise SystemExit(
                f"lake solution is on terrain {lake_terrain[:16]}, mesh is "
                f"{mesh.terrain_hash[:16]}; re-run surface_water.py")
        f_lake = land_fraction_of_class(mesh, grid_export(config), lake)
        depth = (args.lake_dwmax_m if args.lake_dwmax_m is not None
                 else float(model.get("lake_dwmax_m", 0.2)))
        before = float(field[land].mean())
        field = np.where(land, (1.0 - f_lake) * field + f_lake * depth, field)
        lake_report = {
            "source": str(args.lakes),
            "lake_dwmax_m": depth,
            "mean_lake_fraction_of_land_cells": float(f_lake[land].mean()),
            "land_mean_dwmax_before_m": round(before, 5),
            "land_mean_dwmax_after_m": round(float(field[land].mean()), 5),
            "ceiling": "annual lake evaporation is still capped by local "
                       "precipitation; river water never re-enters dwatc",
        }

    output = args.output or (INPUTS / resolution.lower()
                             / f"orogen_{resolution}_surf_{SOIL_WATER_CODE:04d}.sra")
    output.parent.mkdir(parents=True, exist_ok=True)
    write_sra(output, SOIL_WATER_CODE, field)

    weights = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, nlon))
    land_mean = float(np.average(field[land], weights=weights[land])) if land.any() else 0.0

    report = {
        "climatology": rel(args.climatology),
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "code": SOIL_WATER_CODE,
        "field": "dwmax, maximum soil water capacity, metres",
        "lakes": lake_report,
        "soil_map": rel(args.soil_map),
        "soil_map_sha256": hashlib.sha256(args.soil_map.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(args.config.read_bytes()).hexdigest(),
        "resolution": resolution,
        "land_cells": int(land.sum()),
        "land_cells_matched": matched,
        "land_cells_without_soil_left_at_default": unmatched_land,
        "land_mean_capacity_m": land_mean,
        "exoplasim_uniform_default_m": EXOPLASIM_DEFAULT_WSMAX_M,
        "why_this_matters": (
            "landmod.f90 computes runoff as the overflow of this bucket, so a "
            "uniform default makes runoff a namelist property rather than a "
            "soil property. The 0.96 climatology's 2.8% land runoff ratio "
            "against Earth's ~35% is the symptom."),
        "not_enabled_by_default": (
            "Requires model.soil_water_source: pedology in config/planet.yaml. "
            "That key is absent on purpose: adding it moves config_sha256 and "
            "blocks resumption of runs in flight."),
        "output": rel(output),
        "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=PROJECT_ROOT).stdout.strip() or None,
    }
    report_path = output.with_name(output.stem + "_provenance.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n")

    print(f"land cells        {int(land.sum())}, matched {matched}, "
          f"left at default {unmatched_land}")
    print(f"land-mean dwmax   {land_mean:.3f} m against ExoPlaSim's uniform "
          f"{EXOPLASIM_DEFAULT_WSMAX_M} m")
    print(f"\nwrote {rel(output)}")
    print(f"      {report_path.name}")
    print("\nNOT enabled. Set model.soil_water_source: pedology in "
          "config/planet.yaml when re-baselining;\nthat key is absent by design "
          "because adding it moves config_sha256 and blocks resuming runs.")


if __name__ == "__main__":
    main()
