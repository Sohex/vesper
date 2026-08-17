#!/usr/bin/env python3
"""Tectonic and magmatic ore prospectivity, as a 0-1 field per deposit type.

    python minerals/scripts/build_prospectivity.py

Writes `minerals/data/<source_build>/prospectivity.nc` on the native mesh, plus a
provenance report beside it.

## What this is, and what it deliberately is not

It is a FIELD, not deposits. A porphyry system is one to two kilometres against a
15.19 km mesh cell and a vein is under a hundredth of one, so an individual body
is invisible at every resolution this pipeline runs at. Placing one would invent
detail the grid cannot hold. Discrete deposits belong with the downscaling pass,
which is also where glacial overdeepening goes, for the same reason and at almost
the same scale ratio. See `notes/economic-minerals.md`.

It is NOT a lithology and NOT an erodibility modifier. `substrate_class` sets
erodibility and therefore terrain, albedo and therefore climate, texture and
therefore the biosphere, and solutes and therefore the weathering fluxes.
Anything entering it enters the climate path, and a mesh cell would be painted
with ore properties on the strength of a deposit occupying a fraction of a
percent of it.

## Why this is downstream rather than inside Orogen

`notes/economic-minerals.md` assigns tectonic and magmatic genesis to Orogen,
meaning the PROCESSES that concentrate these deposits are the ones Orogen models.
It does not require the arithmetic to happen there, and every input but one is
already exported: craton weight, fold-belt weight, stress, back-arc distance,
substrate and basement class, cover thickness and erosion delta. The missing one,
`lipV`, has its expression in the `flood_basalt` class.

Computing it here buys two things. Iterating on a rule costs no terrain rebuild,
so the layer stays freely rebuildable as the design says it should. And a value
that does not exist during generation cannot feed back into erosion or albedo, so
the constraint above becomes structural rather than a thing to remember.

## What the rules do

They are first order, they live in `config/prospectivity.yaml` rather than here,
and each names the control it encodes. The axis worth understanding is
exhumation, because it is what a rock map alone cannot tell you: a porphyry forms
1-5 km down and is destroyed by deep erosion, while an orogenic gold system forms
5-15 km down and is revealed by it. The same erosion that removes one exposes the
other.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import netCDF4 as nc
import numpy as np
import yaml

from _paths import CONFIG, DATA, PROJECT_ROOT, PROSPECTIVITY

import builds
from orogen import LAND, Export


def normalise(field: np.ndarray, land: np.ndarray) -> np.ndarray:
    """Scale to 0-1 over land, leaving ocean at zero.

    Normalised by the maximum rather than by a quantile: prospectivity is a
    relative statement within this world, and clipping the top would flatten
    exactly the cells the layer exists to identify.
    """
    out = np.zeros_like(field, dtype=np.float32)
    top = field[land].max() if land.any() else 0.0
    if top > 0:
        out[land] = (field[land] / top).astype(np.float32)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-build", default=None,
                    help="override config's source_build")
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()

    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    if args.source_build:
        config["source_build"] = args.source_build
    rules = yaml.safe_load(PROSPECTIVITY.read_text(encoding="utf-8"))

    mesh = Export(builds.mesh_export(config))
    land = mesh.field("surface_class") == LAND
    substrate = mesh.field("substrate_class").astype(int)
    basement = mesh.field("basement_rock").astype(int)
    thickness = mesh.field("cover_thickness").astype(float)
    erosion = mesh.field("erosionDelta").astype(float)
    craton = mesh.field("r_t_craton").astype(float)
    fold = mesh.field("r_t_foldBelt").astype(float)
    stress = mesh.field("r_stress").astype(float)
    stress = stress / stress.max() if stress.max() > 0 else stress
    area = mesh.field("cell_area").astype(float)

    codes = {c["code"]: c["id"] for c in mesh.manifest["lithology"]["rockClasses"]}
    ex = rules["exhumation"]
    cover_ok = thickness >= ex["cover_preserved_km"]
    deep = erosion >= ex["deep_exhumation_delta"]

    fields, summary = {}, {}
    for key, spec in rules["deposits"].items():
        score = np.zeros(land.shape, dtype=float)

        # Host rock. Substrate is what is at the top of the stack; basement is
        # what a deposit sat in. Both count, because a deposit hosted in the
        # basement is still there when cover survives above it.
        for code, weight in (spec.get("hosts") or {}).items():
            if code not in codes:
                continue
            score = np.maximum(score,
                               weight * ((substrate == codes[code])
                                         | (basement == codes[code])))

        if spec.get("fold_belt_weight"):
            score = score * (1.0 + spec["fold_belt_weight"] * fold)
        if spec.get("stress_weight"):
            score = score * (1.0 + spec["stress_weight"] * stress)

        # Exhumation, the axis that separates the hydrothermal families.
        if spec.get("require_cover"):
            score = score * cover_ok
        if spec.get("favour_deep_exhumation"):
            score = score * (1.0 + deep)

        if spec.get("craton_weight"):
            gate = craton >= spec.get("craton_threshold", 0.0)
            score = spec["craton_weight"] * craton * gate

        score[~land] = 0.0
        fields[key] = normalise(score, land)
        lit = land & (fields[key] > 0)
        summary[key] = {
            "name": spec["name"],
            "land_fraction_nonzero": float(area[lit].sum() / area[land].sum()),
            "land_fraction_above_half": float(
                area[land & (fields[key] > 0.5)].sum() / area[land].sum()),
            "mean_over_land": float(np.average(fields[key][land],
                                               weights=area[land])),
        }

    out = args.output or (DATA / str(config["source_build"]) / "prospectivity.nc")
    out.parent.mkdir(parents=True, exist_ok=True)
    with nc.Dataset(out, "w", format="NETCDF4") as data:
        data.createDimension("region", land.size)
        for key, field in fields.items():
            var = data.createVariable(key, "f4", ("region",), zlib=True)
            var[:] = field
            var.long_name = rules["deposits"][key]["name"]
            var.units = "1"
            var.description = ("relative prospectivity, 0-1 over land; NOT a "
                               "deposit and NOT a lithology")
        data.setncattr("vesper_source_build", str(config["source_build"]))
        data.setncattr("vesper_terrain_hash", mesh.terrain_hash)
        data.setncattr("note", "Tectonic and magmatic ore prospectivity. Read "
                               "only downstream of climate; never an input to "
                               "erodibility, albedo or any climate path.")

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source_build": config.get("source_build"),
        "terrain_hash": mesh.terrain_hash,
        "config_sha256": hashlib.sha256(PROSPECTIVITY.read_bytes()).hexdigest(),
        "generator": "minerals/scripts/build_prospectivity.py",
        "git_commit": subprocess.run(["git", "rev-parse", "HEAD"],
                                     capture_output=True, text=True,
                                     cwd=PROJECT_ROOT).stdout.strip() or None,
        "deposits": summary,
        "note": ("Relative within this world, normalised by the land maximum "
                 "per deposit type. A 1.0 is the most favourable cell here, not "
                 "an absolute grade or tonnage, and cross-type comparison of "
                 "the numbers is meaningless."),
    }
    (out.parent / "prospectivity_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"{'deposit':28}{'% land':>9}{'% >0.5':>9}{'mean':>8}")
    for key, s in summary.items():
        print(f"  {s['name'][:26]:26}{100*s['land_fraction_nonzero']:9.2f}"
              f"{100*s['land_fraction_above_half']:9.2f}{s['mean_over_land']:8.3f}")
    print(f"\nwrote {out.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
