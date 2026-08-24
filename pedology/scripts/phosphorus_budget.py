#!/usr/bin/env python3
"""Relative phosphorus geography and basin-concentration hypotheses.

    python pedology/scripts/phosphorus_budget.py

Writes `pedology/analysis/phosphorus_budget.json`.

## Why this exists

`config/pedogenesis.yaml` carries two grounded per-lithology phosphorus rows from
Hartmann et al. (2014) -- how much P a rock CONTAINS, and how much it RELEASES per
unit weathering flux -- and until now **nothing read either of them**. The claim
that this world has P-starved uplands against P-rich basin floors, and that wind
returns some of it, was prose with grounded data sitting beside it and no
computation in between. The result remains a relative diagnostic; ANUT-1
through ANUT-3 own the missing mass-carrying implementation.

## What it computes

Per mesh region, a relative weathering-release score, then the drainage terminal
and geometric basin-floor concentration associated with that score. This ranks
where the outbound leg may be important; it does not carry an absolute P mass,
water flux, dissolution, retention or root-zone delivery.

The aeolian return leg is a source-composition hypothesis rather than a
deposition model. This script does not read `dust_baseline.nc`; it compares the
P content of deflatable source lithologies with the land mean and therefore
cannot settle the sign or magnitude of P delivery at a destination:

**The mapped dust source on this world is its lowest-phosphorus material.** Earth's
Sahara-to-Amazon transport works because the Bodele Depression is a former lake
bed rich in biogenic, P-bearing diatomite. Here the equivalent surface is
evaporite and playa clastic, which Hartmann's own classes put at the bottom of
the P range. That makes a low-P source-composition case, not a dilution result:
deposition can still add P, and its importance depends on deposited mass and
the fraction that dissolves.

**The source-composition hypothesis carries one caveat large enough to invert
its ranking**, and this script reports both sides rather than choosing. Basin
fill is assigned P by lithology
class, at the value Hartmann gives unconsolidated sediment. But basin fill is not
primary rock -- it is whatever the catchment delivered, concentrated by having
nowhere else to go. If closed-basin fill inherits its catchment's P instead of
its class default, the fill is *enriched* rather than depleted and the sign
of the source-composition comparison flips. Which is right is a question about
this world's sediment, not about Hartmann, so both are computed and neither is
presented as a nutrient-delivery answer.
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
from paths import rel  # noqa: E402

sys.path.insert(0, str(PROJECT_ROOT / "lib"))
import builds  # noqa: E402
from orogen import Export, LAND  # noqa: E402
from surface_classes import cover_mask  # noqa: E402

import yaml  # noqa: E402

PEDOGENESIS = HERE.parent / "config" / "pedogenesis.yaml"
OUT = ANALYSIS / "phosphorus_budget.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=OUT)
    ap.add_argument("--surface-classes", type=Path, default=None,
                    help="pedology/analysis/surface_classes.nc. Given, the "
                         "budget reports the two derived classes that carry "
                         "phosphorus with opposite signs: `diatomite` as the "
                         "deflatable source and `pavement` as a sink. SURF-6")
    args = ap.parse_args()

    cfg = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    ped = yaml.safe_load(PEDOGENESIS.read_text(encoding="utf-8"))
    p_ppm = ped["phosphorus_ppm"]
    p_rel = ped.get("phosphorus_release_relative") or {}
    barren = set(cfg["model"]["barren_rock_classes"])

    # The MESH carrier: surface_class, cell_area and substrate_class are
    # native-mesh fields and Export refuses an export without raw/. The
    # literal here read as "the T42 export" and meant "wherever the mesh
    # is". SPAT-2.
    ex = Export(builds.mesh_export(cfg))
    land = ex.surface_class == LAND
    area = ex.cell_area
    rock = ex.substrate_class
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
    # terminal >= 0 is a preserved basin sink; -1 is the ocean; -2 is not
    # land (build_hydrography raises if any land region keeps it).
    endorheic = land & (terminal >= 0)
    exorheic = land & (terminal == -1)

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
        },
        "aeolian_source": {
            "classes": sorted(barren),
            "fraction_of_land": round(float(area[dust_mask].sum() / land_area), 4),
            "mean_content_ppm": round(wmean(content, dust_mask), 1),
            "enrichment_vs_land_mean": round(
                wmean(content, dust_mask) / mean_content, 3),
            "interpretation":
                "Below 1.0 means the mapped source lithology is POORER in "
                "phosphorus than mean land. It does not determine whether a "
                "destination gains or loses P: that also requires emitted and "
                "deposited mass, source provenance, mineral phase and "
                "dissolution. ANUT-3.",
            "caveat_that_could_invert_this":
                "Basin fill is assigned P by lithology class, at Hartmann's "
                "value for unconsolidated sediment. But fill is not primary "
                "rock: it is whatever the catchment delivered, with nowhere "
                "else to go. If it inherits catchment P instead, its source "
                "composition is ENRICHED. Neither case supplies the deposition "
                "or bioavailability needed to settle the return flux. See "
                "catchment_inheritance below.",
        },
        "by_rock_class": {},
    }

    # ---- the two derived classes that move phosphorus, with opposite signs ---
    #
    # SURF-6. Lithology says what the rock IS; these say what the surface has
    # BECOME under this climate, and phosphorus is the budget where the
    # difference has a sign attached.
    #
    # `diatomite` is the deflatable source. The aeolian_source block above reads
    # BARREN LITHOLOGY and concludes the return leg dilutes rather than
    # fertilises -- explicitly noting that Sahara-to-Amazon works because its
    # source is biogenic diatomite rather than evaporite. This world grows
    # diatomite in its strandlines, so the comparison stops being rhetorical:
    # the class is exactly the material the analogy names.
    #
    # `pavement` is a sink, and the sign is the correction Wells et al. (1995)
    # forced. Dust falls through the clast mosaic and accumulates beneath it as
    # a cumulic Av horizon, so phosphorus arriving on a pavement is retained
    # rather than returned. Read as deflation armour it would have been a
    # source; it is the opposite.
    if args.surface_classes is not None:
        derived = {}
        for name, role in (("diatomite", "source"), ("pavement", "sink")):
            # By NAME, from the file's own flag legend. lib/surface_classes.py
            # says why, and raises on a class this build does not carry rather
            # than quietly selecting nothing.
            sel = cover_mask(args.surface_classes, name)
            if sel.shape != rock.shape:
                raise SystemExit(
                    f"surface_cover has {sel.shape[0]} regions and the export "
                    f"has {rock.shape[0]}; they are indexed by region, so a "
                    "mismatch means they were built from different terrain")
            m = land & sel
            derived[name] = {
                "role": role,
                "fraction_of_land": round(float(area[m].sum() / land_area), 4),
                "mean_content_ppm": round(wmean(content, m), 1) if m.any() else None,
                "enrichment_vs_land_mean": (
                    round(wmean(content, m) / mean_content, 3) if m.any() else None),
            }
        derived["note"] = (
            "The two enter the aeolian leg with OPPOSITE signs and neither is a "
            "lithology: diatomite is deflatable biogenic silica, the material "
            "the Sahara-to-Amazon analogue actually runs on, and pavement is a "
            "non-erodible cover over an accretionary horizon that stores what "
            "falls on it. `aeolian_source` above is the lithological source and "
            "does not include either.")
        result["derived_surface_classes"] = derived
    else:
        result["derived_surface_classes"] = {
            "absent": "no --surface-classes given, so the aeolian leg here is "
                      "lithological only: the deflatable diatomite source and "
                      "the pavement sink are both unread. SURF-6."}

    for code, name in sorted(names.items()):
        m = land & (rock == code)
        if not m.any():
            continue
        result["by_rock_class"][name] = {
            "fraction_of_land": round(float(area[m].sum() / land_area), 4),
            "content_ppm": float(p_ppm.get(name, 0.0)),
            "release_relative": float(p_rel.get(name, 0.0)),
        }

    # ---- fill phosphorus as a delivered FLUX, per basin ---------------------
    #
    # Fill is not primary rock and must not be assigned P by lithology class. It
    # is what the catchment delivered, concentrated by a closed basin
    # evaporating the water and keeping the solute. So the enrichment of a
    # basin's floor is a property of its geometry and its catchment's rock, and
    # it has to ARISE rather than be imposed.
    #
    # The test case is the Bodele Depression, the floor of former Lake
    # Mega-Chad and Earth's single largest dust source. What makes it rich is
    # not that it is a basin: it is that a very large catchment was concentrated
    # onto a small floor for long enough that aquatic productivity laid down
    # biogenic sediment, and that it then dried out so wind could deflate it.
    # Large catchment, small floor, currently dry. Each of those is measured
    # here, and a basin qualifies only if it has all three.
    with Dataset(data / "basins.nc") as ds:
        catchment_km2 = np.asarray(ds["catchment_km2"][:])
        floor_km2 = np.asarray(ds["area_at_spill_km2"][:])

    # cell_area is in km2, verified against 4*pi*R^2 for this radius, so the
    # basin floor stays in km2 too. Mixing them was worth a factor of a million
    # and silently produced a catchment-to-floor ratio of zero for every basin.
    delivered = np.zeros(len(catchment_km2))     # relative P score x area
    catch_area = np.zeros(len(catchment_km2))
    order = np.argsort(terminal[land], kind="stable")
    t_sorted = terminal[land][order]
    f_sorted = (flux[land] * area[land])[order]
    a_sorted = area[land][order]
    edges = np.searchsorted(t_sorted, np.arange(-1, len(catchment_km2) + 1))
    for b in range(len(catchment_km2)):
        lo, hi = edges[b + 1], edges[b + 2]
        if hi > lo:
            delivered[b] = f_sorted[lo:hi].sum()
            catch_area[b] = a_sorted[lo:hi].sum()

    floor = np.maximum(floor_km2, 1e-9)          # km2, as cell_area is
    areal_catch = np.divide(delivered, np.maximum(catch_area, 1.0))
    areal_floor = delivered / floor
    enrich = np.divide(areal_floor, areal_catch,
                       out=np.zeros_like(areal_floor), where=areal_catch > 0)
    ratio = np.divide(catch_area, floor,
                      out=np.zeros_like(floor), where=floor > 0)

    # Currently dry? A basin holding water is not a dust source.
    wet = np.zeros(len(catchment_km2), dtype=bool)
    sw = data / "surface_water.nc"
    if sw.is_file():
        with Dataset(sw) as ds:
            if "lake" in ds.variables:
                lake = np.asarray(ds["lake"][:]).astype(bool)
                for b in range(len(catchment_km2)):
                    lo, hi = edges[b + 1], edges[b + 2]
                    if hi > lo and lake[land][order][lo:hi].any():
                        wet[b] = True

    have = catch_area > 0
    bodele = have & (ratio >= 100.0) & ~wet
    strong = have & (ratio >= 20.0) & ~wet
    result["basin_fill_as_flux"] = {
        "note": "Fill phosphorus concentration score per unit floor area, "
                "relative to the same release score spread over its own "
                "catchment; no absolute mass or time unit is implied. Enrichment "
                "equals the catchment-to-floor area ratio when release is "
                "uniform, so it is geometry, not an assumption.",
        "basins_with_catchment": int(have.sum()),
        "catchment_to_floor_ratio": {
            "median": round(float(np.median(ratio[have])), 2),
            "p90": round(float(np.percentile(ratio[have], 90)), 2),
            "p99": round(float(np.percentile(ratio[have], 99)), 2),
            "max": round(float(ratio[have].max()), 2),
        },
        "largest_catchment_km2": round(float(catch_area[have].max())),
        "bodele_analogues": {
            "criteria": "catchment/floor >= 100 and currently dry",
            "count": int(bodele.sum()),
            "earth_reference": "Lake Mega-Chad concentrated roughly 2.5e6 km2 "
                               "onto a roughly 2.5e4 km2 floor, a ratio near 100",
            "max_enrichment": round(float(enrich[bodele].max()), 1)
                              if bodele.any() else None,
        },
        "strongly_concentrating_and_dry": int(strong.sum()),
        "share_of_delivered_p_in_strong_basins": round(
            float(delivered[strong].sum() / max(delivered[have].sum(), 1e-30)), 4),
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
    f = result["basin_fill_as_flux"]
    r = f["catchment_to_floor_ratio"]
    print(f"\ncatchment/floor ratio  median {r['median']}  p99 {r['p99']}  max {r['max']}")
    bo = f["bodele_analogues"]
    print(f"  Bodele analogues (ratio >= 100, dry): {bo['count']}"
          + (f", max enrichment {bo['max_enrichment']}x" if bo['max_enrichment'] else ""))
    print(f"  strongly concentrating and dry (>= 20): {f['strongly_concentrating_and_dry']}")
    print(f"  share of delivered P in those basins: "
          f"{f['share_of_delivered_p_in_strong_basins']:.1%}")
    print(f"\nwrote {rel(args.output)}")


if __name__ == "__main__":
    main()
