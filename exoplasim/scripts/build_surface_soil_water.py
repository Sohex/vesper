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
    python exoplasim/scripts/build_surface_soil_water.py --grid source/<build>/exoplasim-<rung>

**THE GRID IS THE RUNG AND IT IS TOLD, not inferred from a climatology.** The
three rung-dependent inputs are the grid, the ownership mask -- surface code
0172, which `build_boundary_conditions.py` writes at every exported rung -- and
`land_column_states_<res>.txt`, and none of them is a climate. What a
climatology can still do here is FAIL: pass `--climatology` and its axes and its
`lsm` are required to be the grid and the mask this builder already has.

**This reads and does not derive.** The field is
`pedology/data/<build>/land_column_states_<res>.txt`, written by
`pedology/scripts/land_column_properties.py`. It used to be `soilmap.txt`'s
`awc` column, which is pedology's own declared endmember mixture and a
different capacity; WORLD-OF6N settled which of the two the world has.

**TWO CAPACITY COLUMNS, AND THEY ARE TWO QUANTITIES.** The states file carries
both and this script installs exactly one of them:

    awc_mm         plant-available over the whole physical column, field
                   capacity minus the wilting point, which is what a ROOT can
                   remove. LPJ-GUESS installs this and NOTHING here moves it.
    evaporable_mm  the same integral with the SURFACE LAYER's share cut from air
                   dry instead, because a bare drying surface reaches the water
                   between air dry and the wilting point and a root does not.

`evaporable_mm` is what ExoPlaSim installs as `dwmax`, and it is installed only
where the declared water column actually HAS the contract's surface layer at its
top. The two travel together: the extra water is water the top 0.02 m can give
up, so putting it in a column whose top layer is 0.5 m would hand a seasonal
store water that only a drying skin can reach. Where the surface layer is absent
this installs `awc_mm` and cuts every layer from the wilting point, which is what
it did before PHYS-15. `exoplasim/notes/soil-albedo-moisture.md` argues the
floor; `pedology/config/land_column_properties.yaml` declares it.

The identity that keeps the two from being one number written twice is checked
here per cell and not merely inherited: the three layer capacities this script
cuts must sum to the capacity it installs.

**Expect the effect to be small, and build it anyway.** The uniform bucket was
once the leading suspect for this world's land runoff ratio, and an offline
bucket experiment refuted that: shrinking it 3.75-fold moves runoff by 1.06x,
because the bucket sits at a median 15% of capacity, so overflow comes from the
wettest cells and seasons, which saturate at either depth. The low ratio is the
climate, not the namelist. What this field buys is a soil-borne runoff field
rather than a namelist-borne one, which is the loop, and 229 is cheap once the
soil exists. See `pedology/README.md`.

**This changes climate results**, so it belongs in the run whose climatology the
carve verdict will use, not after it. `model.soil_water_source` in
`config/planet.yaml` is the switch and it says which state it is in; the code
defaults to `uniform` when the key is absent, which is what kept runs in flight
resumable while the key was held back, since moving it is a configuration change
`continue_exoplasim.py` refuses to resume across.
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
from paths import rel  # noqa: E402
from provenance import config_stamp  # noqa: E402
from sra import read_sra, write_sra
from builds import grid_export, land_column_states, resolution_of  # noqa: E402
from gridding import (gaussian_grid, grid_geometry,  # noqa: E402
                      land_fraction_of_class, model_label_cells,
                      require_model_labels)
import rungs  # noqa: E402

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

# ExoPlaSim's own default, landmod.f90 `WSMAX_EARTH`. Only cells with NO
# ground keep it now: a cell the binary 0.5 coastline rounding gave to the
# ocean still has a soil column and SPAT-5's land tile needs its capacity,
# so the population is the states file's and not the ownership mask's.
# Where there is genuinely no land this is still meaningless and still has
# to be a sane number.
EXOPLASIM_DEFAULT_WSMAX_M = 0.5

# The model's OWNERSHIP mask at this rung, written by
# build_boundary_conditions.py from the mesh and the declared
# `geography_land_threshold`. It is what the climatology's `lsm` is a copy of --
# index for index, bit-identical, measured in
# notes/audits/grid-convention-and-runoff.md -- and unlike a climatology it
# exists at every exported rung. That is the whole of why this builder is no
# longer pinned to the one rung a run has reached.
LAND_MASK_CODE = 172

# Coordinate rounding shared with pedology/scripts/build_soil.py.
COORD_DECIMALS = 4


def read_land_column_states(path: Path, column_name: str
                            ) -> dict[tuple[float, float], float]:
    """One capacity column in mm, keyed by rounded coordinates.

    From the contract's emitted states, not from the soil map. The file carries
    several columns that would all parse as a capacity and they are NOT
    interchangeable -- `awc_mm` is LPJ-GUESS's and `evaporable_mm` is this
    model's -- so the caller names which one it wants and the file and the
    column are both named in the error.
    """
    lines = path.read_text().splitlines()
    header = lines[0].split()
    if column_name not in header:
        raise SystemExit(
            f"{path} has no {column_name} column. It is the land column "
            "property contract's emitted states; re-emit it with "
            "pedology/scripts/land_column_properties.py.")
    column = header.index(column_name)
    capacity = {}
    for line in lines[1:]:
        parts = line.split()
        capacity[(round(float(parts[0]), COORD_DECIMALS),
                  round(float(parts[1]), COORD_DECIMALS))] = float(parts[column])
    return capacity


def read_retention_states(path: Path) -> dict[tuple[float, float], tuple[float, float]]:
    """Field capacity and wilting point per cell, volumetric.

    The two volumetric states the per-layer capacities are cut between, read
    from the same states file rather than re-derived from texture: the contract
    owns the closure and this script owns no retention curve. The surface
    layer's cut needs the wilting point EXPLICITLY, because its capacity is the
    one that does not stop there.
    """
    lines = path.read_text().splitlines()
    header = lines[0].split()
    for name in ("theta_fc", "theta_wp"):
        if name not in header:
            raise SystemExit(
                f"{path} has no {name} column. The per-layer capacity split is "
                "cut between the two volumetric states; re-emit the states "
                "with pedology/scripts/land_column_properties.py.")
    ifc, iwp = header.index("theta_fc"), header.index("theta_wp")
    out = {}
    for line in lines[1:]:
        parts = line.split()
        out[(round(float(parts[0]), COORD_DECIMALS),
             round(float(parts[1]), COORD_DECIMALS))] = (float(parts[ifc]),
                                                         float(parts[iwp]))
    return out


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


def layer_depth_weights(thicknesses: list[float], physical_m: float,
                        physical_count: int) -> np.ndarray:
    """How much of each physical increment each model water layer owns, in m.

    Returns an array of shape (water layers, physical increments) whose entries
    are overlap DEPTHS, so each column sums to `physical_m` and the whole array
    sums to the column depth.

    The model's water column is a PARTITION of the contract's physical column
    and never an interpolation of it. Where a model boundary lands on an
    increment boundary that is a grouping, which is what LSHY-3's 0.5 and 1.0 m
    cut is: the contract's own 500 mm upper over five increments and 1000 mm
    lower over ten, the only two-layer cut at which the climate column and the
    ecology column share a boundary.

    A boundary may also land INSIDE an increment, which is what the contract's
    0.02 m surface layer does, and splitting by depth there is still a partition
    rather than an interpolation. The reason is the contract's own
    `materials.mineral.vertical_rule`, which is `uniform`: one texture per cell
    applied to every increment, so the retention states are constant within an
    increment and the usable share is one number for it. Half an increment
    therefore holds exactly half its water, by the contract's own rule and not
    by an assumption made here. If that rule ever becomes a profile, this
    function is where the profile has to be integrated instead.

    What is still refused is a column that does not COVER the contract's
    column, because a split over part of it would not sum to the capacity
    `dwmax` carries.
    """
    weights = np.zeros((len(thicknesses), physical_count), dtype=float)
    top = 0.0
    for index, thickness in enumerate(thicknesses):
        if thickness <= 0.0:
            raise SystemExit(
                f"surface.land_water_column.layer_thickness_m[{index}] = "
                f"{thickness} is not positive. It is a depth in metres.")
        bottom = top + thickness
        for cell in range(physical_count):
            lo, hi = cell * physical_m, (cell + 1) * physical_m
            overlap = min(bottom, hi) - max(top, lo)
            if overlap > 0.0:
                weights[index, cell] = overlap
        top = bottom
    column_m = physical_m * physical_count
    if abs(top - column_m) > 1.0e-9:
        raise SystemExit(
            f"the declared water column reaches {top} m of the {column_m} m "
            f"the contract carries. A split over part of the column would not "
            "sum to the capacity dwmax carries.")
    return weights


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--states", type=Path, default=None,
                        help="the land column property contract's emitted "
                             "per-cell states; defaults to the configured "
                             "build's")
    parser.add_argument("--grid", type=Path, default=None,
                        help="the export whose grid this field is cut on; the "
                             "GRID is the rung, so the output path, the "
                             "filename, the land mask and the states file all "
                             "follow it. Defaults to the configured build's")
    parser.add_argument("--land-mask", type=Path, default=None,
                        help=f"staged surface code {LAND_MASK_CODE}, the "
                             "model's ownership mask at this rung; defaults to "
                             "the one build_boundary_conditions.py wrote beside "
                             "the output")
    parser.add_argument("--climatology", type=Path, default=None,
                        help="OPTIONAL CROSS-CHECK and never a carrier: when "
                             "given, its axes and its `lsm` are required to be "
                             "the grid and the mask this builder already has. "
                             "Nothing is read from it")
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

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    model = config["model"]

    # THE GRID IS THE RUNG, and it is TOLD rather than inferred from a
    # climatology. This builder used to take its grid and its land mask off
    # `bootstrap_climatology` and call `require_configured_grid` on it, which
    # made the one rung a run had reached the only rung the field could be cut
    # at -- alone among the thirteen staged surface fields, and with the model
    # refusing to start above it under `soil_water_source: pedology`. What the
    # climatology supplied was the grid and the ownership mask and nothing else;
    # both are properties of the terrain at a rung, and both are already written
    # at every exported rung. WORLD-QGB6.
    #
    # The DETERMINATION edge the climatology also stood for is not lost with it:
    # the climate enters this chain at `pedology/scripts/build_soil.py`, which
    # names its own stage and refuses a soil/climate pairing it did not declare.
    # Reading a climatology here as well put a grid coupling on that argument.
    #
    # Deliberately from the grid, not from config: see builds.resolution_of.
    grid_dir = (args.grid or grid_export(config)).resolve()
    resolution = resolution_of(grid_dir)
    grid_lat, grid_lon, _grid_manifest = grid_geometry(grid_dir)
    nlat, nlon = int(grid_lat.size), int(grid_lon.size)
    ladder_nlat, ladder_nlon, _ = rungs.geometry(resolution)
    if (nlat, nlon) != (ladder_nlat, ladder_nlon):
        raise SystemExit(
            f"{grid_dir} is named {resolution} and ships a {nlat}x{nlon} grid, "
            f"which is not the ladder's {ladder_nlat}x{ladder_nlon} for that "
            "rung. The export and its name disagree; nothing here can say "
            "which is right.")

    if args.states is None:
        # The states file is per build AND RUNG, and the rung here is the
        # GRID's. `builds.land_column_states` reads it out of a config, so it is
        # asked with the grid's rung substituted rather than with the filename
        # spelled a second time -- and `rungs.model_grid` inside it then refuses
        # a rung and a shape that are not one fact.
        states_config = dict(config, model=dict(
            model, resolution=resolution, latitudes=nlat, longitudes=nlon))
        args.states = land_column_states(states_config)

    if not args.states.is_file():
        # WHY THIS IS NOT "run land_column_properties.py" AT EVERY RUNG. That
        # script derives the states from `soilmap_<res>.txt` and nothing else,
        # so the missing file is really a missing soil map, and `build_soil.py`
        # reads a climatology for temperature, precipitation, runoff,
        # evaporation and elevation -- content, not a grid. A soil at a rung
        # therefore needs a CLIMATE at that rung. Naming only the next script up
        # sent a caller at a rung config is not set to round the same loop: it
        # writes the states at CONFIG's rung, this refusal fires again, and
        # nothing says why. WORLD-512R decided the soil waits for the climate
        # rather than being built on a remapped one; the argument and the order
        # that unblocks it are in exoplasim/notes/route-step-criteria.md.
        configured = str(config["model"]["resolution"]).upper()
        raise SystemExit(
            f"{rel(args.states)} does not exist. It is per build AND per rung, "
            "and it is written by pedology/scripts/land_column_properties.py "
            f"from pedology's soil map at {resolution}, which "
            "pedology/scripts/build_soil.py writes from a CLIMATOLOGY at that "
            f"rung.\nconfig/planet.yaml is at {configured}, and neither "
            "pedology script takes a --grid, so both would build at "
            f"{configured} however this step is invoked."
            + ("" if configured == resolution else
               f" Running them now would write the {configured} states and "
               "leave this refusal in place.")
            + f"\nWhat this needs at {resolution}: a climatology of an arm at "
            f"{resolution} on the current staged surface family, declared in "
            "config/planet.yaml, then build_soil.py and "
            "land_column_properties.py at that rung. See "
            "exoplasim/notes/route-step-criteria.md, 'The order that unblocks "
            "the route', and WORLD-CCX6 for the --grid the pedology scripts "
            "lack.")

    # WHICH CAPACITY THIS INSTALLS, and it is decided by the declared GEOMETRY
    # rather than by a switch of its own. The air-dry floor and the surface
    # layer are one selection: `evaporable_mm` is `awc_mm` plus the water the
    # top 0.02 m can give up below the wilting point, so it belongs in a column
    # that has that layer and nowhere else. A column cut at 0.5 m taking
    # `evaporable_mm` would put a drying skin's water into a seasonal store.
    contract = yaml.safe_load(
        (PROJECT_ROOT / "pedology" / "config" / "land_column_properties.yaml")
        .read_text(encoding="utf-8"))
    column = config.get("surface", {}).get("land_water_column", {}) or {}
    physical_m = float(contract["geometry"]["physical_layer_thickness_m"])
    physical_count = int(contract["geometry"]["physical_layer_count"])
    thicknesses = [float(v) for v in column.get(
        "layer_thickness_m", [physical_m * physical_count])]
    surface = contract["surface_layer"]
    surface_thickness_m = float(surface["thickness_m"])
    air_dry = float(surface["air_dry_water_content"])
    has_surface_layer = abs(thicknesses[0] - surface_thickness_m) <= 1.0e-12
    capacity_column = "evaporable_mm" if has_surface_layer else "awc_mm"
    capacity_mm = read_land_column_states(args.states, capacity_column)
    retention = read_retention_states(args.states)

    # THE LAND MASK IS CODE 172 AT THIS RUNG, which is what every other builder
    # in the staged family reads and what the model itself owns land by. The
    # climatology's `lsm` is a copy of this file index for index -- bit-identical,
    # measured in notes/audits/grid-convention-and-runoff.md -- so nothing is
    # given up by reading the mask directly, and a mask that exists at every
    # exported rung is gained.
    mask_path = args.land_mask or (
        INPUTS / resolution.lower()
        / f"orogen_{resolution}_surf_{LAND_MASK_CODE:04d}.sra")
    if not mask_path.is_file():
        raise SystemExit(
            f"{rel(mask_path)} does not exist. It is surface code "
            f"{LAND_MASK_CODE}, the model's ownership mask at {resolution}; run "
            f"exoplasim/scripts/build_boundary_conditions.py --grid {grid_dir} "
            "first.")
    land = read_sra(mask_path, nlat, nlon) > 0.5

    # A CLIMATOLOGY IS A CHECK HERE AND NEVER A CARRIER. Nothing is read from
    # it: what it can do is fail, by being on another rung, by carrying a
    # longitude axis that is not the labels the model writes, or by disagreeing
    # with code 172 about which cells are land. That is the whole of what this
    # builder used to take from one, restated as an assertion so that a caller
    # who has a climatology to hand gets the agreement checked instead of the
    # rung pinned.
    if args.climatology is not None:
        import netCDF4 as nc
        with nc.Dataset(args.climatology) as data:
            require_model_labels(data["lat"][:], data["lon"][:], nlat, nlon,
                                 what=f"{rel(args.climatology)} against the "
                                      f"{resolution} grid")
            clim_land = np.asarray(data["lsm"][0], dtype=float) > 0.5
        if not np.array_equal(clim_land, land):
            raise SystemExit(
                f"{rel(args.climatology)} and {rel(mask_path)} disagree about "
                f"land on {int((clim_land != land).sum())} of {land.size} "
                "cells. The climatology's `lsm` is the mask the model was "
                "handed, read back, so these are two different builds or two "
                "different thresholds and not two opinions to reconcile.")

    # THE STATES FILE CARRIES THE MODEL'S OWN LABELS, and the labels are not the
    # export's. `pedology/scripts/build_soil.py` writes its rows keyed by the
    # longitude a CLIMATOLOGY carries, which runs 0..360 from column 0; the same
    # column is -177.1875 on the export's axis at T21. The two name the same
    # columns and matching one against the other rotates the field by half the
    # planet, agreeing on 0.5955 of cells and refusing nothing. So the labels are
    # resolved to INDICES through the one door for it, which refuses a label
    # that is not on the model's axis instead of snapping it to the nearest
    # column. CLAUDE.md rule 3; lib/gridding.py:model_label_cells.
    keys = list(capacity_mm)
    key_rows, key_cols = model_label_cells(
        np.array([k[1] for k in keys], dtype=float),
        np.array([k[0] for k in keys], dtype=float),
        nlat, nlon, what=f"{rel(args.states)} against the {resolution} grid")
    cell_key: dict[tuple[int, int], tuple[float, float]] = {}
    for key, j, i in zip(keys, key_rows.tolist(), key_cols.tolist()):
        if (j, i) in cell_key:
            raise SystemExit(
                f"{rel(args.states)} puts two soil columns in grid cell "
                f"({j}, {i}): {cell_key[(j, i)]} and {key}. The states file and "
                f"the {resolution} grid are not the same grid.")
        cell_key[(j, i)] = key

    field = np.full((nlat, nlon), EXOPLASIM_DEFAULT_WSMAX_M)
    matched = 0
    unmatched_land = 0
    matched_outside_ownership = 0
    # THE STATES FILE'S POPULATION IS THE GROUND, NOT THE OWNERSHIP MASK.
    # pedology now writes a column for every cell with a positive land
    # fraction, so a cell the binary 0.5 coastline rounding gave to the ocean
    # still has a real capacity and SPAT-5's land tile needs it. Gating this
    # loop on `land` threw those away and left them at the uniform 0.5 m
    # namelist default, which is a placeholder and not a soil property.
    #
    # Inert for the binary model: landmod reads dwmax only under dls > 0.0, so
    # nothing the model integrates today reads a cell this newly fills. `land`
    # is still what says an UNMATCHED cell is an error rather than an absence.
    for j in range(nlat):
        for i in range(nlon):
            key = cell_key.get((j, i))
            value = None if key is None else capacity_mm.get(key)
            if value is None:
                if land[j, i]:
                    unmatched_land += 1
                continue
            field[j, i] = value / 1000.0   # mm -> m, the units dwmax is in
            matched += 1
            if not land[j, i]:
                matched_outside_ownership += 1

    # A RIGHT ANSWER THE OTHER SIDE ALREADY KNOWS, and the one this builder had
    # no way to fail on. pedology writes a soil column for every cell with a
    # positive land fraction and code 172 owns land at `geography_land_threshold`,
    # so the mask's population is a SUBSET of the states file's and every owned
    # land cell must carry a capacity. A cell left at the namelist default is
    # therefore a states file from another build or another rung, or a mask and a
    # states file that were never on one grid -- and the field written from it is
    # a uniform bucket wearing a soil's provenance. The old form counted these
    # and printed the count.
    if unmatched_land:
        raise SystemExit(
            f"{unmatched_land} of {int(land.sum())} cells that surface code "
            f"{LAND_MASK_CODE} owns as land have no soil column in "
            f"{rel(args.states)}, so they would be written at the uniform "
            f"{EXOPLASIM_DEFAULT_WSMAX_M} m namelist default. pedology's "
            "population is every cell with a positive land fraction and the "
            "mask's is a subset of it, so this is a states file and a mask that "
            "are not the same build at the same rung, not an absence of soil.")

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
        from orogen import Export
        from builds import mesh_export
        mesh = Export(mesh_export(config))
        with nc.Dataset(args.lakes) as lds:
            lake = np.asarray(lds["lake"][:]).astype(bool)
            lake_terrain = getattr(lds, "terrain_hash", None)
        if lake_terrain and lake_terrain != mesh.terrain_hash:
            raise SystemExit(
                f"lake solution is on terrain {lake_terrain[:16]}, mesh is "
                f"{mesh.terrain_hash[:16]}; re-run surface_water.py")
        # ONTO THE GRID THIS FIELD IS BEING CUT ON, not the configured one. The
        # lake solution is a property of the mesh and carries no rung; where it
        # lands is decided by the grid the blend is written on.
        f_lake = land_fraction_of_class(mesh, grid_dir, lake)
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
            "source": rel(args.lakes),
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
    # THE SURFACE LAYER'S SHARE IS CUT FROM AIR DRY, and that is the whole of
    # what PHYS-15 changed here. Every layer below it keeps `theta_fc -
    # theta_wp`, what a root can remove; the surface layer takes `theta_fc -
    # air_dry`, what a bare drying surface can. So the split is no longer a
    # normalisation of the usable shares alone -- the two cuts differ per cell by
    # the wilting point -- and the three capacities are built in millimetres and
    # divided by the capacity actually installed at the end.
    shares, usable_count = read_layer_usable_shares(args.states)
    if usable_count != physical_count:
        raise SystemExit(
            f"{args.states} carries {usable_count} usable-share columns and "
            f"the contract declares {physical_count} "
            "physical increments; re-emit the states.")

    # THE GEOMETRY IS THE THICKNESS LIST, and the count is its length. The split
    # is cut at thicknesses, so that list is the whole of what this needs;
    # `run_exoplasim.py` is where a config whose declared count disagrees with
    # the length of its own per-layer lists is refused, and refusing it twice
    # would put the same rule in two places.
    nwater = len(thicknesses)
    depth_weights = layer_depth_weights(thicknesses, physical_m, physical_count)

    # The geometric split is the one a uniform profile gives: thickness over
    # column depth. It is what a cell with no soil takes, and what every land
    # cell would take if `dsoilwf` were the only route -- so it is also the
    # baseline the spread below is reported against. It is `DSOILWF` by the same
    # arithmetic on the same thicknesses, and on a cell that HAS a soil it is
    # wrong at the surface layer by construction, because it cuts by depth where
    # the surface layer's capacity is cut from air dry.
    geometric = np.array(thicknesses, dtype=float)
    geometric = geometric / geometric.sum()

    split = np.tile(geometric[:, None, None], (1, nlat, nlon))
    split_matched = 0
    worst_identity_mm = 0.0
    # Follows the same population as the capacity above: a cell with a soil
    # column gets its own split, whatever the ownership mask calls it.
    for j in range(nlat):
        for i in range(nlon):
            key = cell_key.get((j, i))
            if key is None:
                continue
            usable = shares.get(key)
            states = retention.get(key)
            installed = capacity_mm.get(key)
            if usable is None or states is None or installed is None:
                continue
            theta_fc, theta_wp = states
            # Metres of usable depth per model water layer, which is the only
            # place the geometry enters. Times the volumetric range that layer's
            # capacity is cut over, times 1000, is millimetres of water.
            depth = depth_weights @ usable
            floor = np.full(nwater, theta_wp, dtype=float)
            if has_surface_layer:
                floor[0] = air_dry
            capacities = 1000.0 * depth * (theta_fc - floor)
            total = float(capacities.sum())
            if total <= 0.0:
                continue
            # A RIGHT ANSWER THE CONTRACT ALREADY KNOWS. The three capacities
            # cut here and the one column installed as `dwmax` are the same
            # integral written two ways, so they must agree cell by cell. This
            # is what would catch the surface layer being cut from the wilting
            # point while `evaporable_mm` is installed, which is silent in every
            # other output: the split would still sum to one.
            worst_identity_mm = max(worst_identity_mm, abs(total - installed))
            split[:, j, i] = capacities / total
            split_matched += 1
    # THE TOLERANCE IS THE STATES FILE'S OWN PRECISION and not a physical
    # allowance. It prints capacities to three decimals and volumetric states to
    # six, so the capacity column and the integral rebuilt from the states beside
    # it differ at the file's rounding: measured at 2e-3 mm at the worst cell of
    # this build. 0.01 mm is an order above that and two and a half orders below
    # what the surface layer's floor is worth, so a wilting-point cut installed
    # against an air-dry capacity fails this by a factor of several hundred.
    identity_tolerance_mm = 0.01
    if split_matched and worst_identity_mm > identity_tolerance_mm:
        raise SystemExit(
            f"the per-layer capacities cut here disagree with the "
            f"{capacity_column} column installed as dwmax by "
            f"{worst_identity_mm:.3e} mm at the worst cell, against a tolerance "
            f"of {identity_tolerance_mm:g}. The split and the capacity are one "
            "integral, so this is a floor cut at the wrong state or a column "
            "read by the wrong name, not a rounding error.")

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
        "physical_increments_per_model_layer": [
            round(float(v), 6) for v in depth_weights.sum(axis=1) / physical_m],
        "geometric_split": [round(float(v), 6) for v in geometric],
        "land_cells_matched": split_matched,
        "surface_layer": {
            "present": has_surface_layer,
            "thickness_m": surface_thickness_m,
            "cut_from": "air_dry" if has_surface_layer else "wilting_point",
            "air_dry_water_content": air_dry,
            "note": "the top water layer's capacity is cut from air dry when it "
                    "IS the contract's surface layer, because a bare drying "
                    "surface reaches water a root cannot. Every layer below it "
                    "keeps field capacity minus the wilting point, which is why "
                    "soilsrwp and soilsrfc still describe them.",
        },
        "capacity_identity_residual_mm": worst_identity_mm,
        "capacity_identity_tolerance_mm": identity_tolerance_mm,
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
                      "contract emits, cut onto the declared water column by "
                      "overlap depth and multiplied by the volumetric range "
                      "each layer's capacity is taken over. The same array the "
                      "capacity column is the integral of, so the split and the "
                      "capacity are one arithmetic rather than two that agree, "
                      "and `capacity_identity_residual_mm` is the check.",
        "lakes": "the lake blend moves dwmax and not the split: the lake "
                 "fraction of a cell has no soil profile, so its share of the "
                 "bucket is cut on the surrounding column's shape. The split "
                 "is a property of the soil column and this is the limit of it."
                 if lake_report is not None else None,
        "output": rel(split_output),
        "output_sha256": hashlib.sha256(split_output.read_bytes()).hexdigest(),
    }

    # The quadrature the model's own global budget is taken over, from the one
    # grid constructor. A cosine of the latitude is a second opinion about a
    # weight this module already owns, and on a Gaussian row it is not the same
    # number.
    weights = gaussian_grid(nlat, nlon).cell_area_fraction()
    land_mean = float(np.average(field[land], weights=weights[land])) if land.any() else 0.0

    report = {
        "grid": rel(grid_dir),
        "land_mask": rel(mask_path),
        "climatology": rel(args.climatology) if args.climatology else None,
        "climatology_role": ("cross-check only: the grid is the export's and "
                             "the land mask is surface code "
                             f"{LAND_MASK_CODE}, so a climatology is never a "
                             "carrier here and is not required to build the "
                             "field. WORLD-QGB6."),
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "code": SOIL_WATER_CODE,
        "field": "dwmax, maximum soil water capacity, metres",
        "lakes": lake_report,
        "capacity_split": split_report,
        "land_column_states": rel(args.states),
        "land_column_states_sha256": hashlib.sha256(
            args.states.read_bytes()).hexdigest(),
        "capacity_column": capacity_column,
        "capacity_source": f"the land column property contract's "
                           f"{capacity_column} column, selected by whether the "
                           "declared water column carries the contract's "
                           "surface layer at its top. This script converts mm "
                           "to m and installs it; it derives no capacity of its "
                           "own. `awc_mm` is LPJ-GUESS's and is not what this "
                           "installs when a surface layer is present.",
        "config_sha256": hashlib.sha256(args.config.read_bytes()).hexdigest(),
        "resolution": resolution,
        "land_cells": int(land.sum()),
        "land_cells_matched": matched,
        "cells_matched_outside_the_ownership_mask": matched_outside_ownership,
        "population": ("every cell with a pedology soil column, which is "
                       "cell land fraction > 0 and NOT the binary lsm; the "
                       "binary model reads dwmax only under dls > 0.0"),
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
    inputs = [args.states, mask_path] + [p for p in (args.climatology, args.lakes)
                                        if p is not None]
    report.update(config_stamp(config,
                               "exoplasim/scripts/build_surface_soil_water.py",
                               inputs=inputs))

    report_path = output.with_name(output.stem + "_provenance.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n")

    print(f"land cells        {int(land.sum())} by the ownership mask; "
          f"matched {matched}, of which {matched_outside_ownership} the "
          f"mask calls ocean; left at default {unmatched_land}")
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
