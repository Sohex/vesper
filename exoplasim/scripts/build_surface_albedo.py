#!/usr/bin/env python3
"""Write background land albedo from World Orogen lithology.

ExoPlaSim presets background albedo to a uniform `albland` of 0.22 and reads
codes 174 (broadband), 175 (<0.75 um) and 176 (>0.75 um) if the files exist. On
this planet the substrate is not uniform at all: area-weighted bare-rock albedo
over land is 0.315, and the spread runs 0.10 for basalt to 0.50 for evaporite.
Using 0.22 everywhere is a planetary albedo error of about 0.04, roughly
-12 W/m2 of absorbed flux before snow, cloud and vegetation feedbacks. That is a
large fraction of the entire 0.85-to-0.95 flux sweep, so it is not a detail.

`rock_albedo` in the export exists for this. Its own description calls it "a
surface boundary condition for the climate stage; vegetation and snow override
it where they exist. Orogen does not apply it to anything."

What this is and is not
-----------------------

This is *substrate* albedo, not land-surface albedo. A vegetated Earth-like
surface sits nearer 0.12 to 0.20, well below bare rock. Writing bare rock gives a
genuinely bare-rock planet, which is the honest starting point when no vegetation
model has run yet, but it will run cold relative to a vegetated world, and that
cold bias propagates into whatever vegetation model consumes the result.

`--mode` picks which risk to take:

  lithology  bare-rock albedo as exported. Physically what the surface is before
             anything grows on it. Correct for a first pass whose purpose is to
             feed a vegetation model, and biased cold.
  vegetated  the opposite physical endmember: everything that could carry
             vegetation set to `--vegetation-albedo`, with evaporite left bare
             because nothing grows on a salt pan. Land mean about 0.223. Use it
             with `lithology` to bracket the answer.
  scaled     the lithology *pattern*, rescaled so the land mean matches
             `--target-mean`. Kept for experiments, but not recommended as a
             production setting: it fixes the mean to a guess at the answer, and
             it puts every cell somewhere physically wrong, making basalt darker
             than any real rock and evaporite far too dark for a salt pan.
  uniform    write nothing and let ExoPlaSim default to 0.22.

The two endmembers are about 15 to 19 W/m2 apart in absorbed flux, against 21
W/m2 for the entire 0.85-to-0.95 stellar sweep that produced a 33 K range. The
land albedo is therefore not a bootstrap detail; it is comparable to the largest
forcing this project has varied on purpose, and with a positive
vegetation-albedo feedback on top it may select between distinct equilibria
rather than shifting one. Run both endmembers before trusting either.

Ocean cells are written with the water class value; ExoPlaSim computes ocean
albedo separately, so the value there is inert, but the array has to be full.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from netCDF4 import Dataset
import numpy as np
import yaml

from _paths import CONFIG, INPUTS, SOURCE
from convert_orogen import write_sra
from gridding import land_fraction_of_class, land_weighted
from orogen import Export, LAND

# 174 broadband, 175 below 0.75 um, 176 above. With NSIMPLEALBEDO=0 the
# radiation uses the two-band pair; 174 is written too so the broadband
# diagnostic agrees rather than silently keeping 0.22.
ALBEDO_CODES = (174, 175, 176)

# 212 is forest fraction. It is written here rather than left at its default
# because it is not independent of the albedo assumption: landmod blends snow
# albedo between forested and unforested endpoints by `dforest`, so leaving it at
# ExoPlaSim's uniform 0.5 asserts a half-forested planet while the background
# albedo asserts bare rock. That inconsistency darkens snow on the cold branch,
# which is exactly where the bistability question is decided.
FOREST_CODE = 212

# Forest fraction implied by each albedo mode.
MODE_FOREST_FRACTION = {
    "lithology": 0.0,   # bare rock carries no canopy
    "vegetated": 0.5,   # Earth-like mixed cover where anything grows
    "scaled": 0.25,     # midpoint, matching that mode's compromise character
}


def _rock_id(mesh_dir: Path, code: str) -> int:
    lit = json.loads((mesh_dir / "manifest.json").read_text(encoding="utf-8"))["lithology"]
    for r in lit["rockClasses"]:
        if r["code"] == code:
            return int(r["id"])
    raise KeyError(f"no rock class {code!r} in {mesh_dir}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--mesh", type=Path, default=SOURCE / "exoplasim-T42",
                    help="export carrying raw/; the mesh is the same for all")
    ap.add_argument("--grid", type=Path, default=None,
                    help="export whose grid to target; defaults to the config resolution")
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--mode", choices=("lithology", "vegetated", "scaled", "uniform"),
                    default=None)
    ap.add_argument("--target-mean", type=float, default=0.20,
                    help="land-mean albedo for --mode scaled")
    ap.add_argument("--vegetation-albedo", type=float, default=0.15,
                    help="albedo of vegetated ground for --mode vegetated")
    ap.add_argument("--forest-fraction", type=float, default=None,
                    help="override the forest fraction implied by --mode")
    args = ap.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    model = config["model"]
    mode = args.mode or model.get("land_albedo_source", "lithology")
    nlat, nlon = int(model["latitudes"]), int(model["longitudes"])

    resolution = str(model["resolution"]).upper()
    grid_dir = args.grid or (SOURCE / f"exoplasim-{resolution}")
    output = args.output or (INPUTS / resolution.lower())

    if mode == "uniform":
        for code in ALBEDO_CODES + (FOREST_CODE,):
            p = output / f"orogen_{resolution}_surf_{code:04d}.sra"
            if p.exists():
                p.unlink()
        print("mode=uniform: removed any albedo SRA; ExoPlaSim will use albland=0.22")
        return

    mesh = Export(args.mesh)
    rock = mesh.surface_rock
    evaporite = _rock_id(args.mesh, "evaporite")
    is_land = mesh.surface_class == LAND

    # Per-region substrate albedo, then integrated over land only. Taking this
    # from the gridded export instead would average open water's 0.06 into every
    # coastal cell.
    region_albedo = mesh.rock_albedo.astype(np.float64).copy()
    if mode == "vegetated":
        region_albedo[is_land & (rock != evaporite)] = args.vegetation_albedo

    fraction, alb_grid, _empty = land_weighted(mesh, grid_dir, region_albedo)
    land_cells = fraction >= float(model["geography_land_threshold"])

    raw_fraction, raw_alb, _ = land_weighted(mesh, grid_dir,
                                             mesh.rock_albedo.astype(np.float64))
    area = mesh.cell_area.astype(np.float64)
    raw_mean = float(np.average(mesh.rock_albedo[is_land], weights=area[is_land]))

    if mode == "scaled":
        alb_grid = np.clip(alb_grid * (args.target_mean / raw_mean), 0.05, 0.80)

    # Ocean cells carry the water value; ExoPlaSim computes ocean albedo itself,
    # so it is inert, but the array has to be full.
    water_albedo = float(next(r["albedo"] for r in json.loads(
        (args.mesh / "manifest.json").read_text(encoding="utf-8"))["lithology"]["rockClasses"]
        if r["code"] == "water"))
    field = np.where(land_cells, alb_grid, water_albedo)

    forest_value = (args.forest_fraction if args.forest_fraction is not None
                    else MODE_FOREST_FRACTION[mode])
    vegetable = land_fraction_of_class(mesh, grid_dir, rock != evaporite)
    forest = np.where(land_cells, vegetable * forest_value, 0.0)

    gw = np.fromfile(grid_dir / "grid" / "gauss_weights.bin", dtype="float64")
    def gmean(a):
        num = (np.where(land_cells, a, 0.0).mean(axis=1) * gw).sum()
        den = (land_cells.mean(axis=1) * gw).sum()
        return float(num / den) if den else 0.0
    final_mean = gmean(field)

    output.mkdir(parents=True, exist_ok=True)
    written = []
    for code in ALBEDO_CODES:
        path = output / f"orogen_{resolution}_surf_{code:04d}.sra"
        write_sra(path, code, field)
        written.append(str(path))
    forest_path = output / f"orogen_{resolution}_surf_{FOREST_CODE:04d}.sra"
    write_sra(forest_path, FOREST_CODE, forest)
    written.append(str(forest_path))

    report = {
        "mode": mode,
        "mesh": str(args.mesh),
        "grid": str(grid_dir),
        "resolution": resolution,
        "terrain_hash": mesh.terrain_hash,
        "codes": list(ALBEDO_CODES) + [FOREST_CODE],
        "forest_fraction_value": forest_value,
        "forest_fraction_land_mean": gmean(forest),
        "forest_note": ("dforest blends snow albedo between forested and "
                        "unforested endpoints, so it has to agree with the "
                        "background albedo assumption. ExoPlaSim's default is a "
                        "uniform 0.5."),
        "land_mean_bare_rock": raw_mean,
        "land_mean_written": final_mean,
        "exoplasim_default_albland": 0.22,
        "land_min": float(field[land_cells].min()),
        "land_max": float(field[land_cells].max()),
        "files": written,
        "caveat": ("Substrate albedo, not land-surface albedo. Vegetation and "
                   "snow are applied by the model on top of this."),
    }
    (output / "albedo_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
