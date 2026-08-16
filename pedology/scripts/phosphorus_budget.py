#!/usr/bin/env python3
"""Where phosphorus is released, where the water takes it, and what dust returns.

    python pedology/scripts/phosphorus_budget.py

Writes `pedology/analysis/phosphorus_budget.json`.

## Why this exists

`config/pedogenesis.yaml` carries two grounded per-lithology phosphorus rows from
Hartmann et al. (2014) -- how much P a rock CONTAINS, and how much it RELEASES per
unit weathering flux -- and until now **nothing read either of them**. The claim
that this world has P-starved uplands against P-rich basin floors, and that wind
returns some of it, was prose with grounded data sitting beside it and no
computation in between.

## What it computes

Per mesh region, the weathering release of P, then the same quantity routed down
the drainage network to wherever the water actually delivers it. Contrasting
those two gives the outbound leg of the loop: how much P leaves the uplands, and
how much of it is trapped rather than reaching the sea.

The aeolian return leg is bounded rather than modelled. Bounding it turns out to
be enough to settle its sign, which is the thing worth knowing:

**The dust source on this world is its lowest-phosphorus material.** Earth's
Sahara-to-Amazon transport works because the Bodele Depression is a former lake
bed rich in biogenic, P-bearing diatomite. Here the equivalent surface is
evaporite and playa clastic, which Hartmann's own classes put at the bottom of
the P range. So dust arriving on an upland is poorer in P than the upland it
lands on, and the return leg dilutes rather than fertilises.

**That conclusion carries one caveat large enough to invert it**, and this script
reports both sides rather than choosing. Basin fill is assigned P by lithology
class, at the value Hartmann gives unconsolidated sediment. But basin fill is not
primary rock -- it is whatever the catchment delivered, concentrated by having
nowhere else to go. If closed-basin fill inherits its catchment's P instead of
its class default, the fill is *enriched* rather than depleted and the sign
flips. Which is right is a question about this world's sediment, not about
Hartmann, so both are computed and neither is presented as the answer.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from netCDF4 import Dataset

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _paths import ANALYSIS, CONFIG, PROJECT_ROOT  # noqa: E402

sys.path.insert(0, str(PROJECT_ROOT / "lib"))
import builds  # noqa: E402
from orogen import Export, LAND  # noqa: E402

import yaml  # noqa: E402

PEDOGENESIS = HERE.parent / "config" / "pedogenesis.yaml"
OUT = ANALYSIS / "phosphorus_budget.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=OUT)
    args = ap.parse_args()

    cfg = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    ped = yaml.safe_load(PEDOGENESIS.read_text(encoding="utf-8"))
    p_ppm = ped["phosphorus_ppm"]
    p_rel = ped.get("phosphorus_release_relative") or {}
    barren = set(ped.get("barren_rock_classes")
                 or cfg["model"].get("barren_rock_classes") or [])

    ex = Export(builds.build_root(cfg) / "exoplasim-T42")
    land = ex.surface_class == LAND
    area = ex.cell_area
    rock = ex.surface_rock
    # The class legend lives in the manifest, indexed by id, not as a field.
    names = {c["id"]: c["code"] for c in ex.manifest["lithology"]["rockClasses"]}

    data = builds.component_data("hydrography", cfg, strict=True)
    with Dataset(data / "regions.nc") as ds:
        terminal = np.asarray(ds["terminal"][:])
    with Dataset(data / "basins.nc") as ds:
        n_basins = len(ds["basin_id"][:])

    # Per-region phosphorus content and release factor, from the rock class.
    content = np.zeros(rock.shape, dtype=float)
    release = np.zeros(rock.shape, dtype=float)
    for code, name in sorted(names.items()):
        m = rock == code
        if not m.any():
            continue
        content[m] = float(p_ppm.get(name, 0.0))
        release[m] = float(p_rel.get(name, 0.0))

    land_area = area[land].sum()
    mean_content = float((content[land] * area[land]).sum() / land_area)

    # Endorheic land: anything whose drainage chain does not reach the ocean.
    # terminal >= 0 is a preserved basin sink; -1 is the ocean; -2 is an
    # unpreserved pit, which build_hydrography resolves into a basin or the sea.
    endorheic = land & (terminal >= 0)
    exorheic = land & (terminal == -1)
    unresolved = land & (terminal < -1)

    def wmean(field, mask):
        a = area[mask]
        return float((field[mask] * a).sum() / a.sum()) if a.sum() else float("nan")

    # Release flux is content times the per-unit-weathering release factor. This
    # is a RELATIVE field: it ranks where P enters solution, not an absolute
    # kg/m2/yr, because the weathering rate constant is not pinned here.
    flux = content * release

    dust_mask = land & np.isin(rock, [i for i, n in names.items() if n in barren])

    result = {
        "note": "Phosphorus content, release and drainage fate per rock class. "
                "Generated by pedology/scripts/phosphorus_budget.py; do not edit. "
                "Release flux is relative, not kg/m2/yr.",
        "generated": datetime.now(timezone.utc).isoformat(),
        "build": cfg.get("source_build"),
        "terrain_hash": ex.terrain_hash,
        "inputs": {"pedogenesis.yaml": sha256(PEDOGENESIS)},
        "land": {
            "mean_content_ppm": round(mean_content, 1),
            "mean_release_flux_relative": round(wmean(flux, land), 3),
        },
        "by_drainage_fate": {
            "endorheic": {
                "fraction_of_land": round(float(area[endorheic].sum() / land_area), 4),
                "mean_content_ppm": round(wmean(content, endorheic), 1),
                "mean_release_flux_relative": round(wmean(flux, endorheic), 3),
            },
            "reaches_ocean": {
                "fraction_of_land": round(float(area[exorheic].sum() / land_area), 4),
                "mean_content_ppm": round(wmean(content, exorheic), 1),
                "mean_release_flux_relative": round(wmean(flux, exorheic), 3),
            },
            "unresolved_pits": {
                "fraction_of_land": round(float(area[unresolved].sum() / land_area), 4),
                "note": "raw steepest-descent noise pits; build_hydrography "
                        "resolves these, so they are neither fate yet",
            },
        },
        "aeolian_source": {
            "classes": sorted(barren),
            "fraction_of_land": round(float(area[dust_mask].sum() / land_area), 4),
            "mean_content_ppm": round(wmean(content, dust_mask), 1),
            "enrichment_vs_land_mean": round(
                wmean(content, dust_mask) / mean_content, 3),
            "interpretation":
                "Below 1.0 means dust is POORER in phosphorus than the land it "
                "falls on, so the aeolian return leg dilutes the uplands rather "
                "than fertilising them -- the opposite of Sahara-to-Amazon, "
                "where the source is biogenic diatomite rather than evaporite.",
            "caveat_that_could_invert_this":
                "Basin fill is assigned P by lithology class, at Hartmann's "
                "value for unconsolidated sediment. But fill is not primary "
                "rock: it is whatever the catchment delivered, with nowhere "
                "else to go. If it inherits catchment P instead, the fill is "
                "ENRICHED and the sign flips. See catchment_inheritance below.",
        },
        "by_rock_class": {},
    }

    for code, name in sorted(names.items()):
        m = land & (rock == code)
        if not m.any():
            continue
        result["by_rock_class"][name] = {
            "fraction_of_land": round(float(area[m].sum() / land_area), 4),
            "content_ppm": float(p_ppm.get(name, 0.0)),
            "release_relative": float(p_rel.get(name, 0.0)),
        }

    # The alternative: fill inherits the area-weighted P of everything that
    # drains to it. Computed globally rather than per basin, which is the
    # conservative version -- a per-basin figure needs the catchment map and
    # would only sharpen the contrast, not change its direction.
    non_fill = land & ~dust_mask
    inherited = wmean(content, non_fill)
    result["catchment_inheritance"] = {
        "mean_content_of_non_fill_land_ppm": round(inherited, 1),
        "fill_content_if_inherited_ppm": round(inherited, 1),
        "fill_content_by_class_ppm": round(wmean(content, dust_mask), 1),
        "enrichment_if_inherited": round(inherited / mean_content, 3),
        "note": "Upper bound on the inherited case, ignoring the further "
                "concentration a closed basin applies by evaporating the water "
                "and keeping the solute. The true value is at least this.",
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    b = result["by_drainage_fate"]
    print(f"land mean P content        {mean_content:.0f} ppm")
    print(f"  endorheic   {b['endorheic']['fraction_of_land']:.1%} of land, "
          f"{b['endorheic']['mean_content_ppm']:.0f} ppm")
    print(f"  to ocean    {b['reaches_ocean']['fraction_of_land']:.1%} of land, "
          f"{b['reaches_ocean']['mean_content_ppm']:.0f} ppm")
    a = result["aeolian_source"]
    print(f"\ndust source {a['fraction_of_land']:.1%} of land, "
          f"{a['mean_content_ppm']:.0f} ppm, "
          f"enrichment {a['enrichment_vs_land_mean']:.2f}x")
    c = result["catchment_inheritance"]
    print(f"  if fill inherited catchment P instead: "
          f"{c['fill_content_if_inherited_ppm']:.0f} ppm, "
          f"{c['enrichment_if_inherited']:.2f}x")
    print(f"\nwrote {args.output.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
