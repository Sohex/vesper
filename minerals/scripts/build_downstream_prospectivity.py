#!/usr/bin/env python3
"""Weathering, drainage and brine ore prospectivity, as a 0-1 field per type.

    python minerals/scripts/build_downstream_prospectivity.py

Writes `minerals/data/<source_build>/downstream_prospectivity.nc` on the native
mesh, plus a report beside it.

## Why this is a second script and a second artifact

`notes/economic-minerals.md` splits deposit types by genesis: a deposit belongs
where the process that concentrates it is MODELLED. The tectonic and magmatic
half needs arcs, fold belts, LIPs and cratons, which Orogen has, and
`build_prospectivity.py` computes it from the export alone. This half needs
climate, runoff and the chemical divide, which Orogen does not have.

That difference is not stylistic. The tectonic field is a pure function of the
terrain and survives a re-run of the baseline; this one reads a climatology, a
lake solution and a per-basin brine solve, so it carries a CLIMATE in its
identity as well as a terrain and is regenerated when either moves. Keeping them
in one file would have given the durable half the disposable half's lifetime.

Supergene copper is the case that shows the split is a split and not a handover:
it needs an Orogen porphyry AND a downstream climate, and neither component can
produce it alone. So it multiplies the tectonic file's own field rather than
re-deriving a porphyry from the rock map.

## The constraints, unchanged

It is a FIELD, not deposits. It is NOT a lithology and NOT an erodibility
modifier. Nothing upstream of climate may consume it. `build_prospectivity.py`
and `notes/economic-minerals.md` argue both at length and the argument is the
same here.

## How well grounded each rule is

Unequally, and the config says so per rule rather than leaving a reader to guess.
Two of the brine entries are DERIVED -- soda and gypsum fall out of the chemical
divide with no free choice, because Hardie and Eugster's first calcite decides
irreversibly which way a brine goes and `brine_paths.py` already solves it. The
rest are DECLARED: no bauxite, laterite or placer paper is held in
`references/`, and lithium and boron are not among Meybeck's eight species at
all, so those two rest on a geological association rather than on anything this
pipeline models. The report carries the grounding beside every number.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import netCDF4 as nc
import numpy as np
import yaml

from _paths import CONFIG, DATA, DOWNSTREAM, PROJECT_ROOT

import builds
from gridding import climatology_cells
from orogen import LAND, Export
from paths import climatology_path, rel
from provenance import require_build

EARTH_YEAR_DAYS = 365.2422


def load_module(component: str, name: str):
    """Import another component's script for a function, deliberately.

    Every component ships a private `_paths`, so a plain cross-component import
    hands the imported module OUR ANALYSIS and DATA directories. The swap is
    explicit and is reversed in a `finally`.

    Two functions are borrowed this way and neither is worth a second copy.
    `build_soil.weathering_intensity` is the Walker-Hays-Kasting form normalised
    so Earth's land mean is 1, which is what makes a threshold on it mean
    "times Earth's average". `surface_water.river_discharge` is the leaf-peeling
    accumulation down the drainage tree, which is what turns an upstream source
    rock into a downstream placer.
    """
    root = PROJECT_ROOT / component / "scripts"
    saved = sys.modules.pop("_paths", None)
    sys.path.insert(0, str(root))
    try:
        spec = importlib.util.spec_from_file_location(f"_{component}_{name}",
                                                      root / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(root))
        sys.modules.pop("_paths", None)
        if saved is not None:
            sys.modules["_paths"] = saved
    return module


def normalise(field: np.ndarray, land: np.ndarray) -> np.ndarray:
    """Scale to 0-1 over land, leaving ocean at zero.

    By the maximum rather than by a quantile, exactly as the tectonic file does,
    and for the same reason: prospectivity is a relative statement within this
    world, and clipping the top would flatten the cells the layer exists to find.
    """
    out = np.zeros_like(field, dtype=np.float32)
    top = field[land].max() if land.any() else 0.0
    if top > 0:
        out[land] = (field[land] / top).astype(np.float32)
    return out


def host_weight(spec: dict, substrate, basement, codes) -> np.ndarray:
    """Best host weight per region, from substrate or basement.

    Both count, because a deposit hosted in the basement is still there when
    cover survives above it. Same rule as the tectonic file.
    """
    score = np.zeros(substrate.shape, dtype=float)
    for code, weight in (spec.get("hosts") or {}).items():
        if code not in codes:
            continue
        score = np.maximum(score, weight * ((substrate == codes[code])
                                            | (basement == codes[code])))
    return score


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--climatology", type=Path, default=None)
    ap.add_argument("--brine", type=Path, default=None)
    ap.add_argument("--tectonic", type=Path, default=None,
                    help="prospectivity.nc from build_prospectivity.py; "
                         "supergene copper multiplies its porphyry field")
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()

    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    rules = yaml.safe_load(DOWNSTREAM.read_text(encoding="utf-8"))
    build = str(config["source_build"])

    climatology = args.climatology or climatology_path(root=PROJECT_ROOT)
    if not climatology.is_file():
        raise SystemExit(f"{climatology} does not exist")
    require_build(climatology, "climatology", config)

    hydro = builds.component_data("hydrography", config, strict=True)
    tectonic = args.tectonic or (DATA / build / "prospectivity.nc")
    if not tectonic.is_file():
        raise SystemExit(
            f"{rel(tectonic)} does not exist; run build_prospectivity.py first. "
            "Supergene enrichment needs a porphyry protore and there is no "
            "fallback: a climate window on its own would place copper wherever "
            "the weather suited, which is the error the genesis split exists "
            "to prevent.")
    brine_path = args.brine or (PROJECT_ROOT / "pedology" / "analysis"
                                / "brine_paths.json")
    if not brine_path.is_file():
        raise SystemExit(
            f"{rel(brine_path)} does not exist; run "
            "pedology/scripts/brine_paths.py first. It solves the chemical "
            "divide per basin and that is what decides which evaporite a basin "
            "grows; guessing it would decide the answer.")
    brine = json.loads(brine_path.read_text(encoding="utf-8"))
    if brine.get("source_build") != build:
        raise SystemExit(
            f"{rel(brine_path)} was computed on build "
            f"{brine.get('source_build')!r} and this is {build!r}")

    mesh = Export(builds.mesh_export(config))
    grid_dir = builds.grid_export(config)
    land = mesh.field("surface_class") == LAND
    area = mesh.field("cell_area").astype(float)
    substrate = mesh.field("substrate_class").astype(int)
    basement = mesh.field("basement_rock").astype(int)
    codes = {c["code"]: int(c["id"])
             for c in mesh.manifest["lithology"]["rockClasses"]}

    # --- climate on the mesh ---------------------------------------------
    with nc.Dataset(climatology) as ds:
        clim_lat = np.asarray(ds["lat"][:], dtype=float)
        scale = 1000.0 * 86400.0 * EARTH_YEAR_DAYS
        pr = np.asarray(ds["pr"][:], dtype=float).mean(axis=0) * scale
        evap = -np.asarray(ds["evap"][:], dtype=float).mean(axis=0) * scale
        tas = np.asarray(ds["tas"][:], dtype=float).mean(axis=0) - 273.15
    # Runoff is P - E and NOT `mrro`, which is river-routed net divergence.
    # Clamped at zero: a catchment delivers zero or more, never less.
    runoff = np.maximum(pr - evap, 0.0)

    row, col = climatology_cells(mesh, grid_dir, clim_lat)
    precip = pr[row, col]
    runoff_region = runoff[row, col]
    temperature = tas[row, col]

    soil = load_module("pedology", "build_soil")
    pedo = yaml.safe_load(
        (PROJECT_ROOT / "pedology" / "config" / "pedogenesis.yaml")
        .read_text(encoding="utf-8"))
    intensity = soil.weathering_intensity(runoff_region, temperature,
                                          pedo["weathering"])

    # --- drainage --------------------------------------------------------
    with nc.Dataset(hydro / "regions.nc") as ds:
        terminal = np.asarray(ds["terminal"][:]).astype(np.int64)
        receiver = np.asarray(ds["receiver"][:]).astype(np.int64)
    with nc.Dataset(hydro / "basins.nc") as ds:
        n_basins = ds.dimensions["basin"].size
        basin_ids = [str(x) for x in np.asarray(ds["basin_id"][:])]
    lakes = hydro / "surface_water.nc"
    require_build(lakes, "lake solution", config)
    with nc.Dataset(lakes) as ds:
        discharge = np.asarray(ds["discharge_m3_s"][:], dtype=float)
        fills_to_spill = np.asarray(ds["basin_fills_to_spill"][:]).astype(bool)

    exorheic = land & (terminal == -1)

    with nc.Dataset(tectonic) as ds:
        if ds.getncattr("vesper_terrain_hash") != mesh.terrain_hash:
            raise SystemExit(
                f"{rel(tectonic)} was built on another terrain")
        porphyry = np.asarray(ds["porphyry_cu_mo"][:], dtype=float)
        orogenic = np.asarray(ds["orogenic_au"][:], dtype=float)

    fields: dict[str, np.ndarray] = {}
    summary: dict[str, dict] = {}

    # --- weathering and drainage types ------------------------------------
    for key, spec in rules["deposits"].items():
        if key == "supergene_cu":
            window = ((precip >= spec["minimum_precipitation_mm_yr"])
                      & (precip <= spec["maximum_precipitation_mm_yr"]))
            score = porphyry * window
        elif key == "placer_au":
            water = load_module("hydrography", "surface_water")
            # Gold-bearing catchment area accumulated down the tree, over total
            # catchment area accumulated the same way: the gold-bearing SHARE of
            # everything draining into this reach.
            gold_area = water.river_discharge(mesh, receiver, orogenic * area)
            all_area = water.river_discharge(mesh, receiver,
                                             np.ones_like(area) * area)
            with np.errstate(invalid="ignore", divide="ignore"):
                share = np.where(all_area > 0, gold_area / all_area, 0.0)
            live = land & (discharge > 0)
            cut = (float(np.percentile(discharge[live],
                                       spec["discharge_percentile"]))
                   if live.any() else np.inf)
            score = share * (discharge >= cut)
        else:
            score = host_weight(spec, substrate, basement, codes)
            score = score * (intensity >= spec["minimum_weathering_intensity"])
            score = score * (precip >= spec["minimum_precipitation_mm_yr"])
            if spec.get("require_exorheic"):
                score = score * exorheic
            # Scale by intensity itself, not only by clearing its gate: the
            # difference between three times Earth's mean and ten is the
            # difference between a lateritic soil and an ore body.
            score = score * intensity

        score = np.where(land, score, 0.0)
        fields[key] = normalise(score, land)
        summary[key] = {"name": spec["name"], "grounding": spec["grounding"]}

    # --- brines, per basin then painted onto its floor ---------------------
    index = {b: i for i, b in enumerate(basin_ids)}
    path = np.full(n_basins, "", dtype=object)
    k_conc = np.full(n_basins, np.nan)
    ca_ratio = np.full(n_basins, np.nan)
    for rowj in brine["basins"]:
        i = index[rowj["basin_id"]]
        path[i] = rowj["path"]
        k_conc[i] = rowj["k_ueq_l"]
        ca_ratio[i] = rowj["ca_over_hco3"]
    concentrates = ~fills_to_spill

    # HOW FAR down its path a brine runs, not merely which side it is on. The
    # divide decides the direction irreversibly, but the species in EXCESS at
    # calcite saturation is what dominates every step after, so the size of that
    # excess is the strength of the prediction. A basin at Ca/HCO3 = 0.31 is a
    # soda lake in a way one at 0.98 is not, and collapsing both to "alkaline"
    # would emit a mask wearing the shape of a field.
    with np.errstate(invalid="ignore"):
        alkaline_excess = np.clip(1.0 - ca_ratio, 0.0, None)
        ca_excess = np.clip(ca_ratio - 1.0, 0.0, None)
    alkaline_excess = np.nan_to_num(alkaline_excess, nan=0.0)
    ca_excess = np.nan_to_num(ca_excess, nan=0.0)

    safe = np.clip(terminal, 0, None)
    in_basin = land & (terminal >= 0)
    finite_k = k_conc[np.isfinite(k_conc)]
    k_cut = (float(np.percentile(finite_k,
                                 rules["brines"]["brine_potash"]
                                 ["minimum_k_percentile"]))
             if finite_k.size else np.inf)

    # Catchment host share, for the two entries that rest on an association
    # rather than on the divide. Area-weighted over each basin's catchment.
    def catchment_share(hosts: dict[str, float]) -> np.ndarray:
        weight = np.zeros(substrate.shape, dtype=float)
        for code, w in hosts.items():
            if code in codes:
                weight = np.maximum(weight, w * (substrate == codes[code]))
        valid = in_basin
        tot = np.bincount(terminal[valid], weights=area[valid],
                          minlength=n_basins)
        hit = np.bincount(terminal[valid], weights=(area * weight)[valid],
                          minlength=n_basins)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(tot > 0, hit / np.maximum(tot, 1e-30), 0.0)

    for key, spec in rules["brines"].items():
        per_basin = np.zeros(n_basins)
        if spec.get("path") == "alkaline":
            per_basin = alkaline_excess
        elif spec.get("path") == "ca_rich":
            per_basin = ca_excess
        elif key == "brine_potash":
            resolved = np.asarray(path != "", dtype=bool)
            # Gated on the K percentile, then scaled by the K concentration
            # itself: type IV survival to the end of the sequence is a matter of
            # how much there was to start with.
            per_basin = (resolved
                         & (np.nan_to_num(k_conc, nan=-1.0) >= k_cut)
                         ) * np.nan_to_num(k_conc, nan=0.0)
        elif spec.get("catchment_hosts"):
            per_basin = catchment_share(spec["catchment_hosts"])
            per_basin = np.where(np.asarray(path != "", dtype=bool),
                                 per_basin, 0.0)
        if spec.get("require_no_overflow"):
            per_basin = per_basin * concentrates
        # An evaporite grows on the basin FLOOR, not over its catchment, so the
        # field is painted on the closed-basin fill and nowhere else. Painting
        # the catchment would put salt on the mountains that feed the pan.
        floor = land & mesh.field("is_endorheic").astype(bool) & (terminal >= 0)
        score = np.where(floor, per_basin[safe], 0.0)
        fields[key] = normalise(score, land)
        summary[key] = {"name": spec["name"], "grounding": spec["grounding"]}

    land_area = float(area[land].sum())
    for key, field in fields.items():
        lit = land & (field > 0)
        summary[key].update({
            "land_fraction_nonzero": float(area[lit].sum() / land_area),
            "land_fraction_above_half": float(
                area[land & (field > 0.5)].sum() / land_area),
            "mean_over_land": float(np.average(field[land], weights=area[land])),
        })

    # A rule that never fires is a bug and so is one that fires everywhere. The
    # same audit the derived-surface classifier runs, for the same reason.
    empty = sorted(k for k, s in summary.items()
                   if s["land_fraction_nonzero"] <= 0.0)
    saturated = sorted(k for k, s in summary.items()
                       if s["land_fraction_nonzero"] >= 0.99)

    out = args.output or (DATA / build / "downstream_prospectivity.nc")
    out.parent.mkdir(parents=True, exist_ok=True)
    with nc.Dataset(out, "w", format="NETCDF4") as ds:
        ds.createDimension("region", land.size)
        for key, field in fields.items():
            var = ds.createVariable(key, "f4", ("region",), zlib=True)
            var[:] = field
            var.long_name = summary[key]["name"]
            var.units = "1"
            var.grounding = summary[key]["grounding"]
            var.description = ("relative prospectivity, 0-1 over land; NOT a "
                               "deposit and NOT a lithology")
        ds.setncattr("vesper_source_build", build)
        ds.setncattr("vesper_terrain_hash", mesh.terrain_hash)
        ds.setncattr("climatology", rel(climatology))
        ds.setncattr("note", "Weathering, drainage and brine ore prospectivity. "
                             "Read only downstream of climate; never an input "
                             "to erodibility, albedo or any climate path. "
                             "Carries a CLIMATE in its identity as well as a "
                             "terrain, so it is regenerated when either moves.")

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": "minerals/scripts/build_downstream_prospectivity.py",
        "source_build": build,
        "terrain_hash": mesh.terrain_hash,
        "climatology": rel(climatology),
        "climatology_sha256": hashlib.sha256(climatology.read_bytes()).hexdigest(),
        "tectonic_field": rel(tectonic),
        "brine_paths": rel(brine_path),
        "config_sha256": hashlib.sha256(DOWNSTREAM.read_bytes()).hexdigest(),
        "git_commit": subprocess.run(["git", "rev-parse", "HEAD"],
                                     capture_output=True, text=True,
                                     cwd=PROJECT_ROOT).stdout.strip() or None,
        "deposits": summary,
        "rule_audit": {"never_fires": empty, "always_fires": saturated},
        "note": ("Relative within this world, normalised by the land maximum "
                 "per deposit type. A 1.0 is the most favourable cell here, not "
                 "an absolute grade or tonnage, and cross-type comparison of "
                 "the numbers is meaningless. Read `grounding` before quoting "
                 "any of them: two of the brine entries are derived from the "
                 "chemical divide and the rest are declared judgment."),
        "disposable": ("Driven by a pre-carve climatology, lake solution and "
                       "brine solve. The carve replaces all three. The tectonic "
                       "half in prospectivity.nc does not move with a climate "
                       "and is a separate artifact for that reason."),
    }
    (out.parent / "downstream_prospectivity_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"{'deposit':30}{'grounding':>10}{'% land':>9}{'% >0.5':>9}{'mean':>8}")
    for key, s in summary.items():
        print(f"  {s['name'][:28]:28}{s['grounding']:>10}"
              f"{100*s['land_fraction_nonzero']:9.2f}"
              f"{100*s['land_fraction_above_half']:9.2f}{s['mean_over_land']:8.4f}")
    if empty:
        print(f"\nrules that never fire: {', '.join(empty)}")
    if saturated:
        print(f"rules that fire everywhere: {', '.join(saturated)}")
    print(f"\nwrote {rel(out)}")


if __name__ == "__main__":
    main()
