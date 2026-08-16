#!/usr/bin/env python3
"""Predict each closed basin's brine path from the lithology it drains.

    python pedology/scripts/brine_paths.py

Meybeck (1987) table 2C gives what each rock type RELEASES to solution. In an
exorheic catchment those ions leave, and table 5 tabulates the shares that do.
Here mostly they cannot leave, so a second step decides what happens instead --
and it is the step that decides which evaporite a basin grows.

**The chemical divide** (Hardie and Eugster 1970; Eugster and Jones 1979;
Deocampo and Jones 2014). Evaporative concentration reaches calcite saturation
first. Calcite removes Ca and CO3 in EQUAL EQUIVALENTS, so whichever is in excess
at that moment dominates every later step and the deficient one is driven toward
zero. The first calcite therefore decides, irreversibly, whether the brine
becomes carbonate-rich or carbonate-poor:

    Ca > (HCO3 + CO3)   ->  carbonate-poor: gypsum, then Ca-Cl brine
    (HCO3 + CO3) > Ca   ->  alkaline: Na-CO3 brine, trona and natron

Deocampo and Jones make it quantitative with the Spencer Triangle in
Ca-SO4-(HCO3+CO3), which is three of table 2C's own columns, so the release table
and the fate model compose with nothing in between.

This is a LITHOLOGY-ORDERING result, not a concentration prediction. See the
caveats in the output and at the foot of this file.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import netCDF4 as nc

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from orogen import Export, LAND  # noqa: E402
import builds  # noqa: E402
from paths import rel  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
MEYBECK = ROOT / "pedology" / "data" / "reference" / "meybeck1987_tables.json"

# Orogen rock class -> Meybeck table 2C lithology.
#
# Orogen has 20 classes against Meybeck's 10, so this is a judgment map and it is
# written out in full rather than derived, because a silent default here would
# put every unmatched class on one side of the divide. Each entry says why.
ROCK_TO_MEYBECK = {
    "morb":                "volcanic_rocks",
    "oib":                 "volcanic_rocks",
    "flood_basalt":        "volcanic_rocks",
    "arc_basalt":          "volcanic_rocks",
    "arc_andesite":        "volcanic_rocks",
    "rift_bimodal":        "volcanic_rocks",
    "granite":             "granite",
    "granodiorite":        "granite",
    "gneiss":              "gneiss",
    # Meybeck's class is "gneiss and mica-schists"; schist belongs with it.
    "schist":              "gneiss",
    # Quartzite is metamorphic by origin but chemically a quartz sand: nearly
    # inert. Meybeck's misc_metamorphic is serpentinite/marble/amphibolite and
    # releases 1375 ueq/l Ca, which quartzite emphatically does not.
    "quartzite":           "sandstone",
    "melange":             "misc_metamorphic",
    # "Shelf sandstone / shale" is a mixture. Shale is far more reactive and
    # dominates the solute load of any such mixture, so it is mapped there;
    # --clastic-as-sandstone tests the other choice.
    "shelf_clastic":       "shale",
    "carbonate":           "sedimentary_carbonate_rocks",
    "foreland_clastic":    "shale",
    "continental_clastic": "sandstone",
    "pelagic":             "shale",
    # Orogen's evaporite class is a SALT crust, so halite rather than gypsum.
    # There is no gypsum lithology in this world's rock table at all, which is
    # itself a result: see the note in the output.
    "evaporite":           "halite_evaporite",
    # Playa mud is detrital fill, not an evaporite. It is the host that
    # evaporites grow in, and mapping it to an evaporite would beg the question.
    "playa_clastic":       "shale",
}

SPECIES = ["sio2_umol_l", "ca_ueq_l", "mg_ueq_l", "na_ueq_l", "k_ueq_l",
           "cl_ueq_l", "so4_ueq_l", "hco3_ueq_l"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clastic-as-sandstone", action="store_true",
                    help="map shelf_clastic to sandstone instead of shale, to "
                         "test how much the result leans on that judgment")
    ap.add_argument("--output", type=Path,
                    default=ROOT / "pedology" / "analysis" / "brine_paths.json")
    args = ap.parse_args()

    cfg = builds._config()
    mesh = builds.mesh_export(cfg)
    export = Export(mesh)
    data = builds.component_data("hydrography", cfg, strict=True)

    tables = json.loads(MEYBECK.read_text(encoding="utf-8"))
    t2c = tables["table_2c"]
    cols = {c: i for i, c in enumerate(t2c["columns"])}
    release_table = t2c["rows"]

    mapping = dict(ROCK_TO_MEYBECK)
    if args.clastic_as_sandstone:
        mapping["shelf_clastic"] = "sandstone"

    manifest = json.loads((mesh / "manifest.json").read_text(encoding="utf-8"))
    classes = {c["id"]: c["code"] for c in manifest["lithology"]["rockClasses"]}
    unmapped = sorted({c for i, c in classes.items()
                       if c != "water" and c not in mapping})
    if unmapped:
        raise SystemExit(f"rock classes with no Meybeck mapping: {unmapped}")

    # Per-region release vector, from its surface rock class.
    rock = export.surface_rock.astype(np.int32)
    area = export.cell_area.astype(np.float64)
    land = export.surface_class == LAND
    release = np.zeros((len(SPECIES), export.n_regions))
    for cid, code in classes.items():
        if code == "water":
            continue
        sel = rock == cid
        if not sel.any():
            continue
        row = release_table[mapping[code]]
        for si, sp in enumerate(SPECIES):
            release[si, sel] = row[cols[sp]]

    with nc.Dataset(data / "regions.nc") as d:
        terminal = np.array(d["terminal"][:])
    with nc.Dataset(data / "basins.nc") as d:
        nbasin = d.dimensions["basin"].size
        basin_id = [str(x) for x in np.array(d["basin_id"][:])]
        catchment_km2 = np.array(d["catchment_km2"][:])

    # Area-weighted mean release over each basin's catchment. Weighting by area
    # rather than by runoff is the first-pass simplification; see caveats.
    valid = land & (terminal >= 0)
    w = np.where(valid, area, 0.0)
    tot = np.bincount(terminal[valid], weights=w[valid], minlength=nbasin)
    mean = np.zeros((len(SPECIES), nbasin))
    for si in range(len(SPECIES)):
        s = np.bincount(terminal[valid], weights=(w * release[si])[valid],
                        minlength=nbasin)
        mean[si] = np.where(tot > 0, s / np.maximum(tot, 1e-30), np.nan)

    ca = mean[SPECIES.index("ca_ueq_l")]
    hco3 = mean[SPECIES.index("hco3_ueq_l")]
    so4 = mean[SPECIES.index("so4_ueq_l")]
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = ca / hco3
    resolved = tot > 0
    alkaline = resolved & (ratio < 1.0)
    ca_rich = resolved & (ratio >= 1.0)

    # Spencer Triangle coordinates, normalised to unit sum.
    tri_sum = ca + so4 + hco3
    tri = {"ca": ca / tri_sum, "so4": so4 / tri_sum, "hco3": hco3 / tri_sum}

    area_alk = float(tot[alkaline].sum())
    area_ca = float(tot[ca_rich].sum())
    print(f"basins with a resolved catchment: {int(resolved.sum())} of {nbasin}")
    print(f"  alkaline  (Na-CO3, trona/natron): {int(alkaline.sum()):5}  "
          f"{100*area_alk/(area_alk+area_ca):5.2f}% of catchment area")
    print(f"  Ca-rich   (gypsum / Ca-Cl)      : {int(ca_rich.sum()):5}  "
          f"{100*area_ca/(area_alk+area_ca):5.2f}% of catchment area")
    r = ratio[resolved]
    print(f"  Ca/HCO3 across basins: min {np.nanmin(r):.3f}  median "
          f"{np.nanmedian(r):.3f}  max {np.nanmax(r):.3f}")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "generator": "pedology/scripts/brine_paths.py",
        "source_build": cfg.get("source_build"),
        "terrain_hash": builds.terrain_hash(cfg),
        "meybeck_tables_sha256": hashlib.sha256(
            MEYBECK.read_bytes()).hexdigest(),
        "method": {
            "divide": "calcite removes Ca and CO3 in equal equivalents, so the "
                      "species in excess at calcite saturation dominates all "
                      "subsequent evolution (Hardie and Eugster 1970)",
            "release": "Meybeck (1987) table 2C, per lithology",
            "weighting": "catchment area, NOT runoff",
            "rock_class_map": mapping,
        },
        "caveats": [
            "Area weighting, not runoff weighting. Concentration times runoff is "
            "the flux that actually matters, and runoff is not uniform across a "
            "catchment. This fixes the ORDERING of basins, not their absolute "
            "chemistry. Redo with discharge once a current climatology exists.",
            "Table 2C is temperate-stream release on Earth. The values are not "
            "this world's; they carry the relative behaviour of rock types.",
            "The divide is stated on total equivalents and ignores Mg, which "
            "complicates it through Mg-calcite and dolomite.",
            "Orogen's rock table has NO gypsum lithology. Its evaporite class is "
            "a salt crust, mapped to halite. So a Ca-rich path here can only "
            "arise from halite evaporite or from carbonate-poor silicate mixes, "
            "and predicted gypcrete is correspondingly rare.",
            "Every basin is treated as if it evaporates to saturation. Whether "
            "it actually does is a water balance and belongs to hydrography.",
        ],
        "summary": {
            "basins_total": int(nbasin),
            "basins_resolved": int(resolved.sum()),
            "alkaline_basins": int(alkaline.sum()),
            "ca_rich_basins": int(ca_rich.sum()),
            "alkaline_catchment_fraction": area_alk / (area_alk + area_ca),
            "ca_ratio_min": float(np.nanmin(r)),
            "ca_ratio_median": float(np.nanmedian(r)),
            "ca_ratio_max": float(np.nanmax(r)),
        },
        "basins": [
            {
                "basin_id": basin_id[i],
                "catchment_km2": float(catchment_km2[i]),
                "ca_over_hco3": float(ratio[i]),
                "path": "ca_rich" if ratio[i] >= 1.0 else "alkaline",
                "spencer": {k: float(v[i]) for k, v in tri.items()},
                **{sp: float(mean[si, i]) for si, sp in enumerate(SPECIES)},
            }
            for i in range(nbasin) if resolved[i]
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {rel(args.output)}")


if __name__ == "__main__":
    main()
