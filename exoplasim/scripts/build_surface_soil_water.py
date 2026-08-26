"""Feed the pedology component's soil water capacity back to ExoPlaSim.

This closes the loop the other way. ExoPlaSim's land hydrology is a single
bucket per cell of depth `dwmax`, surface code 0229, and its runoff is literally
that bucket overflowing:

    landmod.f90:1098   drunoff(:) = AMAX1(0., dwatc(:) - dwmax(:)) / deltsec

`configure()` clears every surface field when handed a landmap, so `dwmax` falls
back to the uniform namelist default of `wsmax = 0.5 m` everywhere. That default
is why `model.uniform_land_surface` exists.

The capacity comes from the LAND COLUMN PROPERTY CONTRACT, which derives it
from texture and regolith depth and is the one place in this pipeline that
derives it at all. Supplying it makes runoff a property of the soil rather than
of a namelist, and the soil is in turn a product of the climate and the
biosphere. That is the third loop in this pipeline.

    python exoplasim/scripts/build_surface_soil_water.py

**This reads and does not derive.** The field is
`pedology/data/<build>/land_column_states_<res>.txt`, written by
`pedology/scripts/land_column_properties.py`, whose `awc_mm` column is the
plant-available capacity of the declared physical column: the contract's
retention states at this world's gravity, integrated over the declared column
with the weathered-bedrock rule applied. The same file gives LPJ-GUESS its
per-layer states, so the two land columns agree on how much water there is by
construction. It used to be `soilmap.txt`'s `awc` column, which is pedology's
own declared endmember mixture and a different capacity; WORLD-OF6N settled
which of the two the world has.

**Expect the effect to be small, and build it anyway.** The uniform bucket was
once the leading suspect for this world's land runoff ratio, and an offline
bucket experiment refuted that: shrinking it 3.75-fold moves runoff by 1.06x,
because the bucket sits at a median 15% of capacity, so overflow comes from the
wettest cells and seasons, which saturate at either depth. The low ratio is the
climate, not the namelist. What this field buys is a soil-borne runoff field
rather than a namelist-borne one, which is the loop, and 229 is cheap once the
soil exists. See `pedology/README.md`.

**This changes climate results**, so it belongs in the run whose climatology the
carve verdict will use, not after it. `model.soil_water_source: pedology` is set
in `config/planet.yaml`; the code defaults to `uniform` when the key is absent,
which is what kept runs in flight resumable while the key was held back, since
adding it is a configuration change `continue_exoplasim.py` refuses to resume
across.
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

from _paths import CONFIG, INPUTS, PROJECT_ROOT  # noqa: E402  (puts lib/ on sys.path)
from paths import climatology_path, rel, require_configured_grid  # noqa: E402
from provenance import config_stamp  # noqa: E402
from sra import write_sra

SOIL_WATER_CODE = 229
# The capacity SPLIT of that bucket over the land column's water column, one
# record per layer. WORLD-VJBZ. `dsoilwf` in `landmod_nl` is one shape for the
# whole simulated planet and the split is not one shape: it is the same integral
# `awc_mm` is, cut at the model's layer boundaries instead of summed, so it
# varies with regolith depth and weathered-bedrock fraction exactly as the
# capacity does. Over this build's land cells the upper share of the declared
# 0.5/1.0 m cut runs from the geometric 0.3333 to 0.8882.
#
# INERT AT ONE LAYER, where the split is identically one and the namelist says
# so exactly. It becomes the field the model refuses to run without the moment
# the water column is cut into more than one.
SOIL_WATER_SPLIT_CODE = 2290

# ExoPlaSim's own default, landmod.f90 `WSMAX_EARTH`. Non-land cells keep it,
# since dwmax is meaningless over ocean but must still be a sane number.
EXOPLASIM_DEFAULT_WSMAX_M = 0.5

# Coordinate rounding shared with pedology/scripts/build_soil.py.
COORD_DECIMALS = 4


def read_land_column_states(path: Path) -> dict[tuple[float, float], float]:
    """Plant-available water capacity in mm, keyed by rounded coordinates.

    From the contract's emitted states, not from the soil map. The two carry
    columns that would both parse as a capacity, so the column name is checked
    by name and the file is named in the error.
    """
    lines = path.read_text().splitlines()
    header = lines[0].split()
    if "awc_mm" not in header:
        raise SystemExit(
            f"{path} has no awc_mm column. It is the land column property "
            "contract's emitted states; write it with "
            "pedology/scripts/land_column_properties.py.")
    column = header.index("awc_mm")
    capacity = {}
    for line in lines[1:]:
        parts = line.split()
        capacity[(round(float(parts[0]), COORD_DECIMALS),
                  round(float(parts[1]), COORD_DECIMALS))] = float(parts[column])
    return capacity


def read_layer_usable_shares(path: Path) -> tuple[dict, int]:
    """Per physical layer, the usable share of that layer's capacity, by cell.

    The `u00..` columns of the same states file `awc_mm` comes from. They are
    the weathered-bedrock rule evaluated ONCE, by the contract's own script, and
    the column capacity is this array integrated -- so the split derived from
    them and the capacity installed as `dwmax` are the same arithmetic rather
    than two that agree today.
    """
    lines = path.read_text().splitlines()
    header = lines[0].split()
    columns = [i for i, name in enumerate(header)
               if name.startswith("u") and name[1:].isdigit()]
    if not columns:
        raise SystemExit(
            f"{path} carries no per-layer usable-share columns. The capacity "
            "split is derived from them; re-emit the states with "
            "pedology/scripts/land_column_properties.py.")
    shares = {}
    for line in lines[1:]:
        parts = line.split()
        shares[(round(float(parts[0]), COORD_DECIMALS),
                round(float(parts[1]), COORD_DECIMALS))] = np.array(
                    [float(parts[i]) for i in columns])
    return shares, len(columns)


def layer_grouping(thicknesses: list[float], physical_m: float,
                   physical_count: int) -> list[slice]:
    """Which physical layer set each model water layer is made of.

    The model's water column is a GROUPING of the contract's physical column
    and never a re-cut of it: LSHY-3 chose 0.5 and 1.0 m because that is the
    contract's own 500 mm upper over five of the physical column's own
    increments and 1000 mm lower over ten, which is the only two-layer cut at
    which the climate column and the ecology column share a boundary. This refuses any thickness that is not a
    whole number of those increments, and any set that does not cover the column,
    because either one would make the split an interpolation of a profile rather
    than a partition of it.
    """
    groups = []
    used = 0
    for index, thickness in enumerate(thicknesses):
        count = thickness / physical_m
        if abs(count - round(count)) > 1.0e-9 or round(count) < 1:
            raise SystemExit(
                f"surface.land_water_column.layer_thickness_m[{index}] = "
                f"{thickness} is not a whole number of the contract's "
                f"{physical_m} m physical increment. The model's water column is "
                "a grouping of that column, so a partial layer has no split.")
        groups.append(slice(used, used + round(count)))
        used += round(count)
    if used != physical_count:
        raise SystemExit(
            f"the declared water column covers {used} of the "
            f"{physical_count} increments the contract carries. A split over "
            "part of the column would not sum to the capacity dwmax carries.")
    return groups


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--states", type=Path, default=None,
                        help="the land column property contract's emitted "
                             "per-cell states; defaults to the configured "
                             "build's")
    # Resolved from config.baseline_climatology, not hardcoded. The default
    # here named `climatology_s096` until 2026-08-17: pre-carve terrain under
    # the superseded k2 spectrum. See lib/paths.py:climatology_path.
    parser.add_argument("--climatology", type=Path, default=None,
                        help="supplies the grid and the land mask, so they match "
                             "the soil map exactly; defaults to the configured "
                             "baseline_climatology")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--lakes", type=Path, default=None,
                        help="surface_water.nc; raises dwmax on the lake fraction "
                             "of each cell so those cells reach the wetness "
                             "ceiling and evaporate at the potential rate")
    parser.add_argument("--lake-dwmax-m", type=float, default=None,
                        help="bucket depth on the lake fraction; defaults to "
                             "model.lake_dwmax_m, which is required when "
                             "--lakes is given")
    args = parser.parse_args()
    if args.climatology is None:
        args.climatology = climatology_path()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    model = config["model"]
    if args.states is None:
        import sys as _sys
        _sys.path.insert(0, str(PROJECT_ROOT / "lib"))
        import builds as _b
        args.states = _b.land_column_states(config)
    nlat, nlon = int(model["latitudes"]), int(model["longitudes"])
    resolution = str(model["resolution"]).upper()

    if not args.states.is_file():
        raise SystemExit(
            f"{args.states} does not exist. Run "
            "pedology/scripts/land_column_properties.py.")
    capacity_mm = read_land_column_states(args.states)

    # Grid and land mask come from the climatology, which is the boundary land
    # mask as the model itself saw it, and is the same grid pedology wrote its
    # soil map on. Deriving either independently is how the first version of
    # this script matched zero cells out of 4106.
    import netCDF4 as nc

    with nc.Dataset(args.climatology) as data:
        lat = np.asarray(data["lat"][:], dtype=float)
        lon = np.asarray(data["lon"][:], dtype=float)
        land = np.asarray(data["lsm"][0], dtype=float) > 0.5
    # The same guard `climatology_path` applies, called explicitly because
    # `--climatology` can hand this an arbitrary file that never went through
    # the resolver. One expression, in lib/paths.py.
    require_configured_grid(args.climatology, config)

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
    # exceeds 40% of dwmax (landmod.f90:103-104), so a SHALLOWER bucket saturates
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
    # (landmod.f90:1527), a separate store that discharges to the ocean. So a
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
        # No fallback. The `.get(..., 0.2)` that stood here read as though
        # config supplied the depth, and config had never carried the key at
        # all, so the hardcoded number was what every run got -- confirmed in
        # inputs/t21/orogen_T21_surf_0229_provenance.json.
        if args.lake_dwmax_m is not None:
            depth = float(args.lake_dwmax_m)
        elif "lake_dwmax_m" in model:
            depth = float(model["lake_dwmax_m"])
        else:
            raise SystemExit(
                "config/planet.yaml has no `model.lake_dwmax_m`, and --lakes "
                "needs it: it is the bucket depth written over the lake "
                "fraction of every land cell and therefore part of the field. "
                "Declare it there or pass --lake-dwmax-m.")
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

    # ------------------------------------------------------------------
    # THE CAPACITY SPLIT, code 2290. WORLD-VJBZ.
    #
    # The same integral as `dwmax`, cut at the model's water layer boundaries
    # instead of summed. It is written unconditionally, because at one layer it
    # is identically one and costs a file, and above one layer the model refuses
    # to run without it: the split is a per-cell property and the namelist
    # `dsoilwf` is a shape for the whole simulated planet.
    # ------------------------------------------------------------------
    contract = yaml.safe_load(
        (PROJECT_ROOT / "pedology" / "config" / "land_column_properties.yaml")
        .read_text(encoding="utf-8"))
    physical_m = float(contract["geometry"]["physical_layer_thickness_m"])
    shares, physical_count = read_layer_usable_shares(args.states)
    if physical_count != int(contract["geometry"]["physical_layer_count"]):
        raise SystemExit(
            f"{args.states} carries {physical_count} usable-share columns and "
            f"the contract declares {contract['geometry']['physical_layer_count']} "
            "physical increments; re-emit the states.")

    # THE GEOMETRY IS THE THICKNESS LIST, and the count is its length. The split
    # is cut at thicknesses, so that list is the whole of what this needs;
    # `run_exoplasim.py` is where a config whose declared count disagrees with
    # the length of its own per-layer lists is refused, and refusing it twice
    # would put the same rule in two places.
    column = config.get("surface", {}).get("land_water_column", {}) or {}
    thicknesses = [float(v) for v in column.get("layer_thickness_m",
                                                [physical_m * physical_count])]
    nwater = len(thicknesses)
    groups = layer_grouping(thicknesses, physical_m, physical_count)

    # The geometric split is the one a uniform profile gives: thickness over
    # column depth. It is what a cell with no soil takes, and what every land
    # cell would take if `dsoilwf` were the only route -- so it is also the
    # baseline the spread below is reported against.
    geometric = np.array(thicknesses, dtype=float)
    geometric = geometric / geometric.sum()

    split = np.tile(geometric[:, None, None], (1, nlat, nlon))
    split_matched = 0
    for j in range(nlat):
        for i in range(nlon):
            if not land[j, i]:
                continue
            usable = shares.get((float(lon_signed[i]), float(lat_rounded[j])))
            if usable is None:
                continue
            total = usable.sum()
            if total <= 0.0:
                continue
            split[:, j, i] = [usable[g].sum() / total for g in groups]
            split_matched += 1

    split_output = output.with_name(
        output.name.replace(f"{SOIL_WATER_CODE:04d}",
                            f"{SOIL_WATER_SPLIT_CODE:04d}"))
    write_sra(split_output, SOIL_WATER_SPLIT_CODE, split)

    upper = split[0][land]
    split_report = {
        "code": SOIL_WATER_SPLIT_CODE,
        "field": "dsoilwfc, the capacity of each land water layer as a "
                 "fraction of dwmax, per cell",
        "water_layer_count": nwater,
        "layer_thickness_m": thicknesses,
        "physical_layer_per_model_layer": [g.stop - g.start for g in groups],
        "geometric_split": [round(float(v), 6) for v in geometric],
        "land_cells_matched": split_matched,
        "top_layer_share_over_land": {
            "min": round(float(upper.min()), 4),
            "p50": round(float(np.percentile(upper, 50)), 4),
            "p95": round(float(np.percentile(upper, 95)), 4),
            "max": round(float(upper.max()), 4),
            "mean": round(float(upper.mean()), 4),
            "at_the_geometric_value": round(
                float(np.mean(np.abs(upper - geometric[0]) < 1.0e-6)), 4),
        },
        "derivation": "the per-layer usable shares the land column property "
                      "contract emits, grouped onto the declared water column "
                      "and normalised. The same array the awc_mm column is the "
                      "integral of, so the split and the capacity are one "
                      "arithmetic rather than two that agree.",
        "lakes": "the lake blend moves dwmax and not the split: the lake "
                 "fraction of a cell has no soil profile, so its share of the "
                 "bucket is cut on the surrounding column's shape. The split "
                 "is a property of the soil column and this is the limit of it."
                 if lake_report is not None else None,
        "output": rel(split_output),
        "output_sha256": hashlib.sha256(split_output.read_bytes()).hexdigest(),
    }

    weights = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, nlon))
    land_mean = float(np.average(field[land], weights=weights[land])) if land.any() else 0.0

    report = {
        "climatology": rel(args.climatology),
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "code": SOIL_WATER_CODE,
        "field": "dwmax, maximum soil water capacity, metres",
        "lakes": lake_report,
        "capacity_split": split_report,
        "land_column_states": rel(args.states),
        "land_column_states_sha256": hashlib.sha256(
            args.states.read_bytes()).hexdigest(),
        "capacity_source": "the land column property contract's awc_mm column. "
                           "This script converts mm to m and installs it; it "
                           "derives no capacity of its own.",
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
            "soil property; the model's low land runoff ratio is the "
            "symptom."),
        "enabled_by": (
            "model.soil_water_source in config/planet.yaml (currently set to "
            "pedology). Changing the key is a configuration change "
            "continue_exoplasim.py refuses to resume across, so moving it "
            "blocks runs in flight."),
        "output": rel(output),
        "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=PROJECT_ROOT).stdout.strip() or None,
    }
    # The fourth staged field, and it had the WEAK half of the stamp: a hash of
    # config/planet.yaml, which cannot separate an edited comment from an edited
    # parameter, and which nothing was checking. `config_stamp` keeps the parsed
    # config beside it so `config_drift` has something to compare.
    # `inputs` covers what `source_config` cannot: the derived files this
    # generator read. The bespoke `*_sha256` keys above recorded some of the
    # same hashes and nothing compared them; `check_consistency.py` compares
    # these. lib/provenance.py:input_stamp carries the argument.
    inputs = [args.states] + [p for p in (args.climatology, args.lakes)
                              if p is not None]
    report.update(config_stamp(config,
                               "exoplasim/scripts/build_surface_soil_water.py",
                               inputs=inputs))

    report_path = output.with_name(output.stem + "_provenance.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n")

    print(f"land cells        {int(land.sum())}, matched {matched}, "
          f"left at default {unmatched_land}")
    print(f"land-mean dwmax   {land_mean:.3f} m against ExoPlaSim's uniform "
          f"{EXOPLASIM_DEFAULT_WSMAX_M} m")
    upper_stats = split_report["top_layer_share_over_land"]
    print(f"capacity split    {nwater} layer(s); top share over land "
          f"{upper_stats['min']:.4f} to {upper_stats['max']:.4f}, p50 "
          f"{upper_stats['p50']:.4f} against the geometric "
          f"{geometric[0]:.4f}")
    print(f"\nwrote {rel(output)}")
    print(f"      {rel(split_output)}")
    print(f"      {report_path.name}")
    # Say what is actually true rather than printing the same warning forever.
    # This read "NOT enabled" unconditionally, long after the key was set, which
    # tells a reader to go and do something already done.
    if str(config["model"].get("soil_water_source", "uniform")) != "uniform":
        print(f"\nENABLED: model.soil_water_source is "
              f"{config['model']['soil_water_source']!r}, so code 229 reaches "
              "the model and replaces ExoPlaSim's uniform 0.5 m field capacity.")
    else:
        print("\nNOT enabled. Set model.soil_water_source: pedology in "
              "config/planet.yaml when re-baselining; adding that key moves "
              "config_sha256 and blocks resuming existing runs.")


if __name__ == "__main__":
    main()
