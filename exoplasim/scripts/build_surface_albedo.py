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
             because nothing grows on a salt pan. Land mean about 0.197. Use it
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

# 174 broadband, 175 below 0.75 um, 176 above. With NSIMPLEALBEDO=0 the
# radiation uses the two-band pair; 174 is written too so the broadband
# diagnostic agrees rather than silently keeping 0.22.
ALBEDO_CODES = (174, 175, 176)


def _rock_id(export: Path, code: str) -> int:
    lit = json.loads((export / "manifest.json").read_text(encoding="utf-8"))["lithology"]
    for r in lit["rockClasses"]:
        if r["code"] == code:
            return int(r["id"])
    raise KeyError(f"no rock class {code!r} in {export}")


def _surface_rock(export: Path) -> np.ndarray:
    with Dataset(export / "planet.nc") as ds:
        return np.asarray(ds["surface_rock"][:])


def load_albedo(export: Path, nlat: int, nlon: int):
    with Dataset(export / "planet.nc") as ds:
        alb = np.asarray(ds["rock_albedo"][:], dtype=np.float64)
        sc = np.asarray(ds["surface_class"][:])
        area = np.asarray(ds["grid_cell_area"][:], dtype=np.float64)
    if alb.shape != (nlat, nlon):
        raise RuntimeError(
            f"{export} is {alb.shape}, config asks for {(nlat, nlon)}; "
            "use the export whose grid matches the model resolution"
        )
    return alb, sc, area


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--export", type=Path, default=SOURCE / "exoplasim-T42")
    ap.add_argument("--output", type=Path, default=INPUTS / "t42")
    ap.add_argument("--mode", choices=("lithology", "vegetated", "scaled", "uniform"),
                    default=None)
    ap.add_argument("--target-mean", type=float, default=0.20,
                    help="land-mean albedo for --mode scaled")
    ap.add_argument("--vegetation-albedo", type=float, default=0.15,
                    help="albedo of vegetated ground for --mode vegetated")
    args = ap.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    model = config["model"]
    mode = args.mode or model.get("land_albedo_source", "lithology")
    nlat, nlon = int(model["latitudes"]), int(model["longitudes"])

    if mode == "uniform":
        for code in ALBEDO_CODES:
            p = args.output / f"orogen_T42_surf_{code:04d}.sra"
            if p.exists():
                p.unlink()
        print("mode=uniform: removed any albedo SRA; ExoPlaSim will use albland=0.22")
        return

    alb, sc, area = load_albedo(args.export, nlat, nlon)
    land = sc == 1
    raw_mean = float(np.average(alb[land], weights=area[land]))

    field = alb.copy()
    if mode == "vegetated":
        # Everything that can carry vegetation does. Evaporite is left at its
        # bare value: a playa stays a playa, and it is 20.2% of this planet's
        # land, which puts a hard floor under how dark the world can get.
        rock = _surface_rock(args.export)
        evaporite = _rock_id(args.export, "evaporite")
        field[land & (rock != evaporite)] = args.vegetation_albedo
    elif mode == "scaled":
        # Scale about the land mean so the pattern is preserved and the mean
        # moves to the target. Clipped to a physical range afterwards.
        field[land] = np.clip(alb[land] * (args.target_mean / raw_mean), 0.05, 0.80)

    final_mean = float(np.average(field[land], weights=area[land]))
    args.output.mkdir(parents=True, exist_ok=True)
    written = []
    for code in ALBEDO_CODES:
        path = args.output / f"orogen_T42_surf_{code:04d}.sra"
        write_sra(path, code, field)
        written.append(str(path))

    report = {
        "mode": mode,
        "export": str(args.export),
        "terrain_hash": json.loads(
            (args.export / "manifest.json").read_text(encoding="utf-8")
        )["hashes"]["finalElevation"],
        "codes": list(ALBEDO_CODES),
        "land_mean_bare_rock": raw_mean,
        "land_mean_written": final_mean,
        "exoplasim_default_albland": 0.22,
        "land_min": float(field[land].min()),
        "land_max": float(field[land].max()),
        "files": written,
        "caveat": ("Substrate albedo, not land-surface albedo. Vegetation and "
                   "snow are applied by the model on top of this."),
    }
    (args.output / "albedo_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
