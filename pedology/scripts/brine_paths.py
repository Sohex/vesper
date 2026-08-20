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

**Weighted by discharge, not by area.** What reaches a basin is concentration
times water, and runoff on this world spans four orders of magnitude inside a
single large catchment, so the two weightings are not a refinement apart. A dry
interior contributes area and almost no solute; a wet rim contributes the
chemistry. `--weighting area` retains the first pass for comparison, and the
report carries how far apart the two land.

Weight is LOCAL runoff generation, not the river discharge passing a region.
Discharge is the accumulated sum down the tree, so weighting by it would count
every upstream region again at every downstream one and let a handful of cells
near the sink decide a basin. `solute_routing.py` says so where it computes it.

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
from provenance import require_build  # noqa: E402

from _paths import climatology_path  # noqa: E402
import solute_routing  # noqa: E402
from solute_routing import SPECIES  # noqa: E402

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
    # Melange -> SHALE, corrected 2026-08-16, and the correction mattered more
    # than any other entry in this map.
    #
    # It was `misc_metamorphic`, chosen when the melange rule could not fire and
    # the class was 0.00% of land, so nothing rested on it. The arc fix made the
    # rule reachable and melange is now over 6% of land, at which point the
    # choice was supplying close to half of this planet's silicate CO2 drawdown.
    #
    # It was wrong twice over. Meybeck's misc_metamorphic is MARBLE: 2.3% of
    # Earth's outcrop, Ca 1375 ueq/l against sedimentary carbonate's 2560, and
    # Ca + Mg = 1740 against HCO3 1730 -- a near-exact carbonate balance, which
    # is the signature of carbonate dissolution rather than silicate weathering.
    # So it gave a forearc province marble chemistry, AND that bicarbonate was
    # entering the silicate total at full weight when most of it is
    # rock-derived.
    #
    # What Orogen's class actually is decides the replacement. It is not the
    # accretionary prism alone: `elevation.js` assigns it to the whole forearc
    # province -- prism, forearc basin and serpentinite together -- because at
    # ~15 km cells they cannot be separated. That province is greywacke,
    # argillite and clastic basin fill by volume, which is Meybeck's shale, and
    # is what shelf, foreland and pelagic clastics already map to.
    #
    # `gneiss` is not the alternative it looks like. A forearc is not felsic
    # crystalline basement, and mapping it there would swap one wrong rock for
    # another while happening to give a smaller number.
    #
    # KNOWN UNDER-REPRESENTATION: serpentinite. Melanges carry it and Meybeck
    # has a peridotite row for exactly that rock. It is left out because the
    # serpentinite fraction of a forearc is volumetrically minor, nothing here
    # constrains it, and the effect is small: at a generous 10% of melange the
    # planet's silicate CO2 total moves about half a percent.
    #
    # The bias directions, since intuition gets one of them backwards.
    # Meybeck's peridotite is LOWER in bicarbonate than shale, 450 against 580,
    # so omitting serpentinite biases CO2 and Ca HIGH, not low. It biases Mg
    # (500 against 240) and silica (180 against 150) LOW. Ultramafic rock being
    # an efficient CO2 sink per unit area is a statement about reaction rate,
    # not about the solute concentration this table carries.
    "melange":             "shale",
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

def single_lithology_check(export, terminal, mapping, release_table, cols,
                           n_basins, mean, weight):
    """Basins draining exactly one rock class must return that rock's own row.

    A definitional identity, and the one thing here that can FAIL. The weighted
    mean of a constant is that constant whatever the weights are, so a basin
    whose whole catchment is one lithology has a known right answer that does not
    depend on the weighting at all. Any indexing slip between region, cell and
    basin, any misalignment of the climate join, any mix-up of species columns
    breaks it; a wrong-but-plausible weighting does not survive it.

    Comparing the two weightings against each other could only ever tell us they
    differ. `docs/src/practice/failure-modes.md` class 17.
    """
    rock = export.substrate_class.astype(np.int32)
    classes = {c["id"]: c["code"]
               for c in export.manifest["lithology"]["rockClasses"]}
    valid = (terminal >= 0) & (weight > 0)
    order = np.argsort(terminal[valid], kind="stable")
    tb = terminal[valid][order]
    rb = rock[valid][order]
    bounds = np.searchsorted(tb, np.arange(n_basins + 1))

    worst = 0.0
    tested = 0
    worst_basin = None
    for b in range(n_basins):
        lo, hi = bounds[b], bounds[b + 1]
        if hi - lo < 1:
            continue
        first = rb[lo]
        if not (rb[lo:hi] == first).all():
            continue
        tested += 1
        row = release_table[mapping[classes[int(first)]]]
        for si, sp in enumerate(SPECIES):
            expect = row[cols[sp]]
            got = mean[si, b]
            err = abs(got - expect) / max(abs(expect), 1e-12)
            if err > worst:
                worst, worst_basin = err, b
    return {"single_lithology_basins": tested,
            "worst_relative_error": float(worst),
            "worst_basin_index": None if worst_basin is None else int(worst_basin)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--clastic-as-sandstone", action="store_true",
                    help="map shelf_clastic to sandstone instead of shale, to "
                         "test how much the result leans on that judgment")
    ap.add_argument("--weighting", choices=("discharge", "area"),
                    default="discharge",
                    help="which weighting decides the reported path. Both are "
                         "always computed and reported; this picks the primary.")
    ap.add_argument("--climatology", type=Path, default=None,
                    help="climatology the runoff weighting is taken from; "
                         "defaults to config.baseline_climatology")
    ap.add_argument("--output", type=Path,
                    default=ROOT / "pedology" / "analysis" / "brine_paths.json")
    args = ap.parse_args()

    cfg = builds._config()
    mesh = builds.mesh_export(cfg)
    grid_dir = builds.grid_export(cfg)
    export = Export(mesh)
    data = builds.component_data("hydrography", cfg, strict=True)
    # NOT resolved: `rel()` and the recorded path should name the climatology
    # where this repo keeps it, and resolving follows symlinks out of the tree.
    climatology = args.climatology or climatology_path()
    if not climatology.is_file():
        raise SystemExit(f"{climatology} does not exist")
    require_build(climatology, "climatology", cfg)

    tables = json.loads(MEYBECK.read_text(encoding="utf-8"))
    t2c = tables["table_2c"]
    cols = {c: i for i, c in enumerate(t2c["columns"])}
    release_table = t2c["rows"]

    mapping = dict(ROCK_TO_MEYBECK)
    if args.clastic_as_sandstone:
        mapping["shelf_clastic"] = "sandstone"

    release = solute_routing.region_release(export, mapping, release_table, cols)

    area = export.cell_area.astype(np.float64)
    land = export.surface_class == LAND
    runoff_m3, runoff_mm, _, _ = solute_routing.region_runoff_m3_yr(
        export, grid_dir, climatology)

    with nc.Dataset(data / "regions.nc") as d:
        terminal = np.array(d["terminal"][:])
    with nc.Dataset(data / "basins.nc") as d:
        nbasin = d.dimensions["basin"].size
        basin_id = [str(x) for x in np.array(d["basin_id"][:])]
        catchment_km2 = np.array(d["catchment_km2"][:])

    weights = {
        "discharge": np.where(land, runoff_m3, 0.0),
        "area": np.where(land, area, 0.0),
    }
    results = {}
    for name, w in weights.items():
        results[name] = solute_routing.basin_totals(terminal, w, release, nbasin)

    mean, total, tot_w = results[args.weighting]
    area_w = results["area"][2]

    def paths(mean_, tot_):
        ca_ = mean_[SPECIES.index("ca_ueq_l")]
        hco3_ = mean_[SPECIES.index("hco3_ueq_l")]
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio_ = ca_ / hco3_
        return ratio_, tot_ > 0

    ratio, resolved = paths(mean, tot_w)
    ratio_area, resolved_area = paths(results["area"][0], results["area"][2])
    ratio_disch, resolved_disch = paths(results["discharge"][0],
                                        results["discharge"][2])

    so4 = mean[SPECIES.index("so4_ueq_l")]
    ca = mean[SPECIES.index("ca_ueq_l")]
    hco3 = mean[SPECIES.index("hco3_ueq_l")]
    tri_sum = ca + so4 + hco3
    with np.errstate(invalid="ignore", divide="ignore"):
        tri = {"ca": ca / tri_sum, "so4": so4 / tri_sum, "hco3": hco3 / tri_sum}

    alkaline = resolved & (ratio < 1.0)
    ca_rich = resolved & (ratio >= 1.0)
    # Area is the comparable denominator whichever weighting decided the path: a
    # basin's share of the endorheic landscape is its catchment, not its water.
    area_alk = float(area_w[alkaline].sum())
    area_ca = float(area_w[ca_rich].sum())

    # A basin whose whole catchment has P <= E receives no water and therefore no
    # solute, so the divide is UNDETERMINED there rather than defaulting to the
    # area answer. Reporting it as alkaline because the rocks would have been
    # alkaline had it rained is exactly the silent fallback that
    # `docs/src/practice/failure-modes.md` class 2 is about.
    dry = (~resolved_disch) & resolved_area
    flipped = resolved_disch & resolved_area & (
        (ratio_disch >= 1.0) != (ratio_area >= 1.0))

    identity = single_lithology_check(export, terminal, mapping, release_table,
                                      cols, nbasin, mean,
                                      weights[args.weighting])
    identity["criterion"] = "worst relative error below 1e-9"
    identity["passes"] = bool(identity["worst_relative_error"] < 1e-9)

    # Silica delivered per basin: umol/l x m3/yr -> mol/yr. 1 m3 is 1000 l and
    # umol is 1e-6 mol, so the two conversions cancel to 1e-3.
    si = SPECIES.index("sio2_umol_l")
    silica_mol_yr = results["discharge"][1][si] * 1e-3

    r = ratio[resolved]
    print(f"weighting: {args.weighting}  (climatology {rel(climatology)})")
    print(f"basins with a resolved catchment: {int(resolved.sum())} of {nbasin}")
    print(f"  alkaline  (Na-CO3, trona/natron): {int(alkaline.sum()):5}  "
          f"{100*area_alk/(area_alk+area_ca):5.2f}% of catchment area")
    print(f"  Ca-rich   (gypsum / Ca-Cl)      : {int(ca_rich.sum()):5}  "
          f"{100*area_ca/(area_alk+area_ca):5.2f}% of catchment area")
    print(f"  Ca/HCO3 across basins: min {np.nanmin(r):.3f}  median "
          f"{np.nanmedian(r):.3f}  max {np.nanmax(r):.3f}")
    print(f"  basins with a catchment but no runoff: {int(dry.sum())}")
    print(f"  basins whose path flips between weightings: {int(flipped.sum())}")
    print(f"  single-lithology identity: {identity['single_lithology_basins']} "
          f"basins, worst relative error {identity['worst_relative_error']:.2e} "
          f"-> {'PASS' if identity['passes'] else 'FAIL'}")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "generator": "pedology/scripts/brine_paths.py",
        "source_build": cfg.get("source_build"),
        "terrain_hash": builds.terrain_hash(cfg),
        "climatology": rel(climatology),
        "climatology_sha256": hashlib.sha256(climatology.read_bytes()).hexdigest(),
        "meybeck_tables_sha256": hashlib.sha256(MEYBECK.read_bytes()).hexdigest(),
        "method": {
            "divide": "calcite removes Ca and CO3 in equal equivalents, so the "
                      "species in excess at calcite saturation dominates all "
                      "subsequent evolution (Hardie and Eugster 1970)",
            "release": "Meybeck (1987) table 2C, per lithology",
            "weighting": args.weighting,
            "weighting_definition":
                "local runoff generation P - E per mesh region, clamped at zero, "
                "times region area. NOT the accumulated river discharge, which "
                "would count every upstream region again at every downstream one.",
            "climate_join": "gridding.climatology_cells: longitude index for "
                            "index, latitude by nearest Gaussian centre",
            "rock_class_map": mapping,
        },
        "identity_check": identity,
        "weighting_comparison": {
            "note": "The two weightings are two answers, not a test of either. "
                    "What is here is how far apart they land, so a result quoted "
                    "from one is quoted knowing what the other said. The identity "
                    "check above is the part that can fail.",
            "basins_resolved_discharge": int(resolved_disch.sum()),
            "basins_resolved_area": int(resolved_area.sum()),
            "basins_with_catchment_but_no_runoff": int(dry.sum()),
            "catchment_area_without_runoff_km2": float(catchment_km2[dry].sum()),
            "basins_whose_path_flips": int(flipped.sum()),
            "alkaline_catchment_fraction_discharge": float(
                area_w[resolved_disch & (ratio_disch < 1.0)].sum()
                / max(area_w[resolved_disch].sum(), 1e-30)),
            "alkaline_catchment_fraction_area": float(
                area_w[resolved_area & (ratio_area < 1.0)].sum()
                / max(area_w[resolved_area].sum(), 1e-30)),
        },
        "caveats": [
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
            "A basin whose catchment generates no runoff under this climatology "
            "is reported UNDETERMINED rather than falling back to the area "
            "answer. Both weightings are in the per-basin rows either way.",
            "Runoff comes from a pre-carve climatology. The carve removes closed "
            "basins, darkens the land and warms the world, so every number here "
            "is regenerated after the next terrain iteration.",
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
            "endorheic_silica_mol_per_year": float(
                np.nansum(silica_mol_yr[resolved_disch])),
        },
        "basins": [
            {
                "basin_id": basin_id[i],
                "catchment_km2": float(catchment_km2[i]),
                "runoff_m3_per_year": float(results["discharge"][2][i]),
                "silica_mol_per_year": float(silica_mol_yr[i]),
                "ca_over_hco3": float(ratio[i]),
                "ca_over_hco3_discharge": float(ratio_disch[i]),
                "ca_over_hco3_area": float(ratio_area[i]),
                "path": ("undetermined_no_runoff" if not resolved[i]
                         else "ca_rich" if ratio[i] >= 1.0 else "alkaline"),
                "spencer": {k: float(v[i]) for k, v in tri.items()},
                **{sp: float(mean[si_, i]) for si_, sp in enumerate(SPECIES)},
            }
            for i in range(nbasin) if resolved_area[i]
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {rel(args.output)}")


if __name__ == "__main__":
    main()
