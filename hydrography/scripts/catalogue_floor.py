#!/usr/bin/env python3
"""Where the endorheic basin catalogue's floor is, and what sets it.

    python hydrography/scripts/catalogue_floor.py
    python hydrography/scripts/catalogue_floor.py --build precarve-craton-10m

Worldbuilding. Vesper is an invented planet; every quantity here is a property
of a World Orogen export of its terrain, measured on the artifact.

The carve list is drawn from the preserved basin set, and loop A closes by
Orogen consuming that list, so whatever decides which depressions enter the
catalogue decides the terrain of every later pass. `selectBasins` in
`vendor/orogen/js/basins.js` applies three floors and then a nesting rule. This
script asks which of them actually decides the answer on a given build, and
whether the answer is a property of the world or of the mesh it was sampled on.

## What "physical" means here, and why the three floors do not share it

`minAreaKm2` and `minDepthKm` are declared as absolute landform sizes: the same
planet should preserve the same basins at any region count. `minCells` is
declared as a resolution floor, a signal-to-noise test that a depression
spanning a handful of cells is a mesh artifact whatever its area works out to.
So a catalogue is PHYSICS-LIMITED when the absolute floors decide it and
RESOLUTION-LIMITED when the cell count does.

That framing has a third case, and this script measures it. The depth floor is
compared against the depression's depth in the model's DIMENSIONLESS elevation
parameter, not in kilometres, while being named and published in kilometres.
The model-unit-to-km curve is strongly nonlinear above sea level and linear
below it, so one model-unit threshold is many different physical depths
depending on where the depression sits. A floor like that is neither physical
nor resolution-set; it is unit-set, and no region count moves it.

## The controls, which have right answers

A measurement of "which floor binds" is worth nothing unless the reconstruction
it is made from is the selection Orogen actually performed. Three checks, each
able to fail:

  REPRODUCTION  Reapplying the published criteria and the nesting rule to the
                published catalogue must return the published preserved set,
                id for id.

  SENSITIVITY   Loosening every floor must admit strictly more and tightening
                every floor strictly fewer. This is what makes the reproduction
                evidence: a routine that ignored the criteria and echoed the
                published list would pass REPRODUCTION and fail here. A single
                floor is deliberately not asserted on, because a floor another
                floor dominates moves nothing when nudged, and that inertness is
                the measurement rather than a fault.

  RESIDUE       The below-sea-level counts this script derives from the raw
                fields must equal the generator's own `landSeaMask` counts in
                the manifest. Two independent paths to one number.

## The reduction semantics, named

An area is EXTENSIVE and is summed over `cell_area`; a share is that sum over
the land-area sum; a count share is a different statistic and is reported
beside it rather than in place of it. `lib/gridding.py` owns the rule. The
distinction is load-bearing here: refining the mesh multiplies the COUNT of
missed depressions while leaving their total AREA nearly unchanged.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import platform
import subprocess
from pathlib import Path

import numpy as np

from _paths import ANALYSIS, PROJECT_ROOT, SOURCE  # noqa: F401

from builds import mesh_export_of  # noqa: E402  (lib/ on the path by _paths)
from orogen import Export  # noqa: E402

# Builds registered in `lib/orogen.py` that carry a payload here. Named rather
# than globbed so a half-linked worktree fails loudly instead of measuring
# whatever happens to be on disk.
#
# The catalogue is a property of the MESH, and the mesh is stored once per build
# under whichever export directory carries `raw/`. `mesh_export_of` IDENTIFIES
# that carrier instead of naming the rung it happens to live under, which is the
# distinction `lib/builds.py` exists to keep: a rung literal in a path is how a
# consumer ends up reading a grid nothing chose.
DEFAULT_BUILDS = ("precarve-craton", "precarve-craton-10m")


def select(catalogue: list[dict], min_depth: float, min_area: float,
           min_cells: int, depth_key: str = "depth") -> set[int]:
    """`selectBasins`'s threshold path, including the nesting rule.

    A transcription of `vendor/orogen/js/basins.js`, deliberately literal: the
    ordering is by (nestDepth, index) because the nesting filter asks whether an
    ANCESTOR is already chosen and that question only has a stable answer if
    outer basins are considered first. `depth_key` exists so the counterfactual
    below can ask what the same selection would do with the depth floor applied
    to the physical depth instead of the model-unit one.
    """
    parent = [c["parentIndex"] for c in catalogue]
    qualifies = [
        c[depth_key] >= min_depth
        and c["areaKm2"] >= min_area
        and c["cellCount"] >= min_cells
        for c in catalogue
    ]
    order = sorted(range(len(catalogue)),
                   key=lambda i: (catalogue[i]["nestDepth"], i))
    chosen: set[int] = set()
    for i in order:
        if not qualifies[i]:
            continue
        p = parent[i]
        nested = False
        while p >= 0:
            if p in chosen:
                nested = True
                break
            p = parent[p]
        if nested:
            continue
        chosen.add(i)
    return chosen


def controls(catalogue: list[dict], preserved_ids: set[str],
             min_depth: float, min_area: float, min_cells: int) -> dict:
    """Reproduction and sensitivity. Either failing makes everything below void."""
    rebuilt = select(catalogue, min_depth, min_area, min_cells)
    rebuilt_ids = {catalogue[i]["id"] for i in rebuilt}
    reproduces = rebuilt_ids == preserved_ids

    # A reproduction that ignored the criteria altogether would still return the
    # published set, so reproduction alone is not evidence. The sensitivity
    # check moves ALL THREE floors at once, in each direction, and requires the
    # answer to move with them. Both directions, because a reproduction that
    # hardcoded the published list would be insensitive to either.
    #
    # A SINGLE floor is deliberately not asserted on. A floor another floor
    # dominates is genuinely inert on that build, moving nothing when nudged,
    # and that inertness is the measurement rather than a fault -- it is exactly
    # what "which floor decides this catalogue" means. So each floor is nudged
    # alone and the result REPORTED, and only the joint perturbation is a gate.
    looser = select(catalogue, min_depth / 2.0, min_area / 2.0,
                    max(1, min_cells - 1))
    tighter = select(catalogue, min_depth * 2.0, min_area * 2.0, min_cells * 2)
    alone = {
        "min_cells_minus_one": select(catalogue, min_depth, min_area,
                                      max(1, min_cells - 1)),
        "min_area_halved": select(catalogue, min_depth, min_area / 2.0,
                                  min_cells),
        "min_depth_halved": select(catalogue, min_depth / 2.0, min_area,
                                   min_cells),
    }
    return {
        "reproduces_published_selection": reproduces,
        "reproduced_count": len(rebuilt_ids),
        "published_count": len(preserved_ids),
        "only_in_reproduction": len(rebuilt_ids - preserved_ids),
        "only_in_published": len(preserved_ids - rebuilt_ids),
        "looser_admits_more": len(looser) > len(rebuilt),
        "tighter_admits_fewer": len(tighter) < len(rebuilt),
        "sensitive_to_criteria": len(looser) > len(rebuilt) > len(tighter),
        "floor_moves_the_answer_when_relaxed_alone": {
            k: len(v) != len(rebuilt) for k, v in alone.items()},
        "note": "The reproduction is the instrument. It passing means the floor "
                "accounting below is about the selection Orogen performed. "
                "sensitive_to_criteria is what makes that pass evidence: a "
                "reproduction blind to the criteria would return the same set "
                "under both perturbations. A floor that moves nothing when "
                "relaxed alone is inert on this build in the relaxing "
                "direction, which is a finding and not a fault.",
    }


def ladder(catalogue: list[dict], chosen: set[int], area: np.ndarray) -> dict:
    """Detected, qualifying, preserved: where the population goes."""
    return {
        "detected": len(catalogue),
        "preserved": len(chosen),
        "preserved_area_km2": float(area[sorted(chosen)].sum()),
        "largest_preserved_km2": float(area[sorted(chosen)].max()),
        "smallest_preserved_km2": float(area[sorted(chosen)].min()),
        "median_preserved_km2": float(np.median(area[sorted(chosen)])),
    }


def floors(catalogue: list[dict], chosen: set[int], area: np.ndarray,
           min_depth: float, min_area: float, min_cells: int) -> dict:
    """What each floor costs, measured by relaxing it and re-selecting.

    Counting rejections directly would overstate every floor, because the
    nesting rule then drops many of the depressions a relaxed floor admits. The
    honest question is what the PRESERVED set would have been, so each floor is
    relaxed and the whole selection re-run.
    """
    depth = np.array([c["depth"] for c in catalogue])
    cells = np.array([c["cellCount"] for c in catalogue])
    pass_depth, pass_area, pass_cells = (depth >= min_depth,
                                         area >= min_area,
                                         cells >= min_cells)

    def relaxed(**kw) -> dict:
        s = select(catalogue,
                   kw.get("min_depth", min_depth),
                   kw.get("min_area", min_area),
                   kw.get("min_cells", min_cells))
        return {
            "preserved": len(s),
            "added_vs_as_built": len(s) - len(chosen),
            "added_fraction_of_as_built": (len(s) - len(chosen)) / max(len(chosen), 1),
            "preserved_area_km2": float(area[sorted(s)].sum()),
            "added_area_km2": float(area[sorted(s)].sum() - area[sorted(chosen)].sum()),
        }

    return {
        "rejections_before_nesting": {
            "fails_cells_alone": int((pass_depth & pass_area & ~pass_cells).sum()),
            "fails_area_alone": int((pass_depth & pass_cells & ~pass_area).sum()),
            "fails_depth_alone": int((pass_area & pass_cells & ~pass_depth).sum()),
            "passes_all_three": int((pass_depth & pass_area & pass_cells).sum()),
        },
        "relax_min_cells": relaxed(min_cells=1),
        "relax_min_area": relaxed(min_area=0.0),
        "relax_min_depth": relaxed(min_depth=-1.0e9),
        "note": "relax_min_cells is the mesh's remaining bite on this build: "
                "the basins the region count alone is keeping out. It is the "
                "quantity the resolution-limited claim is about.",
    }


def depth_floor_units(catalogue: list[dict], chosen: set[int],
                      min_depth: float, min_area: float,
                      min_cells: int) -> dict:
    """The depth floor is compared in model units and published in kilometres.

    So the physical floor it enforces is not one number. This reports the range
    of physical depths the single declared threshold corresponds to, how much of
    the preserved set is shallower than the declared value, and what the
    selection would be if the same threshold were applied to `depthKm`.
    """
    depth = np.array([c["depth"] for c in catalogue])
    depth_km = np.array([c["depthKm"] for c in catalogue])
    sel = sorted(chosen)

    # Depressions sitting within a hair of the threshold: what physical depth
    # did the declared number actually demand of them?
    at_threshold = np.abs(depth - min_depth) < 0.02 * min_depth
    in_km = select(catalogue, min_depth, min_area, min_cells,
                   depth_key="depthKm")

    return {
        "declared_value": min_depth,
        "declared_units_in_manifest": "km",
        "units_actually_compared": "dimensionless model elevation",
        "physical_depth_at_the_threshold_km": {
            "samples": int(at_threshold.sum()),
            "min": float(depth_km[at_threshold].min()) if at_threshold.any() else None,
            "median": float(np.median(depth_km[at_threshold])) if at_threshold.any() else None,
            "max": float(depth_km[at_threshold].max()) if at_threshold.any() else None,
        },
        "preserved_depth_km": {
            "min": float(depth_km[sel].min()),
            "p5": float(np.percentile(depth_km[sel], 5)),
            "median": float(np.median(depth_km[sel])),
        },
        "preserved_shallower_than_declared": int((depth_km[sel] < min_depth).sum()),
        "preserved_shallower_than_declared_fraction": float((depth_km[sel] < min_depth).mean()),
        "counterfactual_floor_in_km": {
            "preserved": len(in_km),
            "net_vs_as_built": len(in_km) - len(chosen),
            "admitted_that_are_not_preserved_now": len(in_km - chosen),
            "dropped_that_are_preserved_now": len(chosen - in_km),
        },
        "note": "A counterfactual, not a proposal. Changing the comparison "
                "changes the preserved set and therefore the terrain, which is "
                "a loop A decision and costs a generation.",
    }


def below_sea_level_residue(export: Export, manifest: dict) -> dict:
    """Rule 1's number: dry closed-basin floor the catalogue does not carry.

    `surface_class` and `land_mask` disagree exactly over land below sea level.
    Most of that sits inside a preserved basin and `is_endorheic` finds it; the
    rest is depressions that did not clear the selection floors, which is why
    reconstructing land as `land_mask | is_endorheic` loses them. This measures
    how much, BY AREA, because that is the quantity a share of land means.
    """
    surface_class = export.field("surface_class")
    land_mask = export.field("land_mask")
    is_endorheic = export.field("is_endorheic")
    cell_area = export.field("cell_area")
    elevation_km = export.field("elevation_km")

    land = surface_class == 1
    land_area = float(cell_area[land].sum())          # extensive
    disagree = land & (land_mask == 0)
    inside = disagree & (is_endorheic != 0)
    outside = disagree & (is_endorheic == 0)

    published = manifest.get("landSeaMask", {})
    return {
        "land_area_km2": land_area,
        "below_sea_level_dry_land": {
            "regions": int(disagree.sum()),
            "area_km2": float(cell_area[disagree].sum()),
            "area_fraction_of_land": float(cell_area[disagree].sum() / land_area),
        },
        "inside_a_preserved_basin": {
            "regions": int(inside.sum()),
            "area_km2": float(cell_area[inside].sum()),
        },
        "outside_any_preserved_basin": {
            "regions": int(outside.sum()),
            "area_km2": float(cell_area[outside].sum()),
            "area_fraction_of_land": float(cell_area[outside].sum() / land_area),
            "area_share_of_the_disagreement": float(
                cell_area[outside].sum() / cell_area[disagree].sum()),
            "count_share_of_the_disagreement": float(outside.sum() / disagree.sum()),
            "deepest_km": float(elevation_km[outside].min()),
        },
        "endorheic_footprint": {
            "area_km2": float(cell_area[is_endorheic != 0].sum()),
            "area_fraction_of_land": float(
                cell_area[is_endorheic != 0].sum() / land_area),
            "note": "The depression FOOTPRINT, which is not the endorheic "
                    "CATCHMENT. The catchment share is "
                    "manifest.basins.drainageConsistency.fractionOfLand and is "
                    "the larger number by roughly a factor of three.",
        },
        "control_against_manifest": {
            "regions_agree": int(disagree.sum()) == published.get("disagreeingRegions"),
            "inside_agree": int(inside.sum()) == published.get("disagreeingInsidePreservedBasin"),
            "published_regions": published.get("disagreeingRegions"),
            "published_inside": published.get("disagreeingInsidePreservedBasin"),
        },
    }


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def measure(label: str, root: Path) -> dict:
    export = Export(root)
    manifest = export.manifest
    basins = manifest["basins"]
    resolution = basins["resolution"]
    catalogue = basins["catalogue"]
    preserved_ids = {p["id"] for p in basins["preserved"]}

    min_depth = float(resolution["minDepthKm"])
    min_area = float(resolution["minAreaKm2"])
    min_cells = int(resolution["minCells"])
    area = np.array([c["areaKm2"] for c in catalogue])

    checks = controls(catalogue, preserved_ids, min_depth, min_area, min_cells)
    if not checks["reproduces_published_selection"]:
        raise SystemExit(
            f"{label}: the selection could not be reproduced from the published "
            "criteria, so nothing below would be a measurement of it")
    if not checks["sensitive_to_criteria"]:
        raise SystemExit(
            f"{label}: loosening and tightening every floor left the selection "
            "where it was, so the reproduction is not reading the criteria and "
            "reproducing the published set proves nothing")

    residue = below_sea_level_residue(export, manifest)
    agreement = residue["control_against_manifest"]
    if not (agreement["regions_agree"] and agreement["inside_agree"]):
        raise SystemExit(
            f"{label}: the below-sea-level residue derived from the raw fields "
            f"disagrees with the manifest's own landSeaMask counts, so one of "
            "the two is reading a different mask")

    chosen = select(catalogue, min_depth, min_area, min_cells)
    return {
        "label": label,
        "export": str(root.relative_to(PROJECT_ROOT)) if root.is_absolute() else str(root),
        "terrain_hash": export.terrain_hash,
        "catalogue_hash": export.catalogue_hash,
        # `threshold` means no drainage hypothesis was applied, which on this
        # project's lineage means pre-carve. Reported as what the manifest says
        # rather than as a derived boolean, because a preserve list that
        # happened to preserve everything is not the same thing.
        "selection_source": basins.get("selectionSource"),
        "pre_carve_note": "A PRE-CARVE build is a LIMIT, not a state. Every "
                          "number here is what sits in source/ before any carve "
                          "verdict has been applied.",
        "resolution": resolution,
        "controls": checks,
        "ladder": ladder(catalogue, chosen, area),
        "floors": floors(catalogue, chosen, area, min_depth, min_area, min_cells),
        "depth_floor_units": depth_floor_units(catalogue, chosen, min_depth,
                                               min_area, min_cells),
        "below_sea_level_residue": residue,
        "manifest_sha256": _sha256(root / "manifest.json"),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--build", action="append", default=None,
                    help="a registered build to measure; repeatable. Default is "
                         "every build in DEFAULT_BUILDS.")
    ap.add_argument("--output", type=Path,
                    default=ANALYSIS / "catalogue_floor.json")
    args = ap.parse_args()

    names = args.build or list(DEFAULT_BUILDS)
    targets = [(name, mesh_export_of(SOURCE / name)) for name in names]

    builds = []
    for label, root in targets:
        if not (root / "manifest.json").is_file():
            raise SystemExit(f"no manifest at {root}; link the payload first")
        print(f"measuring {label} at {root}")
        builds.append(measure(label, root))

    for b in builds:
        r, f, d = b["resolution"], b["floors"], b["depth_floor_units"]
        res = b["below_sea_level_residue"]["outside_any_preserved_basin"]
        print(f"\n{b['label']}: {r['numRegions']:,} regions, "
              f"mean cell {r['avgCellAreaKm2']:.1f} km2")
        print(f"  detected {b['ladder']['detected']:,}  "
              f"preserved {b['ladder']['preserved']:,}")
        print(f"  relaxing minCells would preserve "
              f"{f['relax_min_cells']['preserved']:,} "
              f"({f['relax_min_cells']['added_fraction_of_as_built']:+.1%})")
        print(f"  preserved shallower than the declared depth floor: "
              f"{d['preserved_shallower_than_declared']:,} "
              f"({d['preserved_shallower_than_declared_fraction']:.1%})")
        print(f"  below-sea-level dry land outside any preserved basin: "
              f"{res['area_km2']:,.0f} km2 "
              f"({res['area_fraction_of_land']:.4%} of land) over "
              f"{res['regions']:,} regions")

    out = {
        "note": "Modelled terrain of an invented planet. A pre-carve build is a "
                "LIMIT, not a state: every basin number here is what sits in "
                "source/ before any carve verdict has been applied.",
        "generated": datetime.datetime.now(datetime.UTC).isoformat(),
        "provenance": {
            "script": "hydrography/scripts/catalogue_floor.py",
            "git_commit": subprocess.run(
                ["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT,
                capture_output=True, text=True).stdout.strip() or None,
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
        "builds": builds,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {args.output}")


if __name__ == "__main__":
    main()
