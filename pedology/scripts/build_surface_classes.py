#!/usr/bin/env python3
"""Derived surface classes: what the surface has BECOME, on two independent axes.

    python pedology/scripts/build_surface_classes.py

World Orogen gives primary lithology: what the rock IS. It does not give what
this world's climate and drainage have made of it, and that second thing is what
albedo, dust emission and the phosphorus cycle key on. It is derivable here
rather than upstream because it depends on precipitation, evaporation and lake
extent that Orogen never sees.

## Two axes, not one chain

`surface_cover` is what wind and light see. `duricrust` is what is cementing at
or below the surface. A cell can be loess over calcrete, and that is a landscape
rather than a contradiction to be resolved.

The `carved-zoned-v2` build lost basin fill because a single lithology cover
chain decided which deposit sat on top by the order the branches happened to be
written, and 203 preserved basins ended up with no fill cell anywhere. The lesson
is not "be careful with branch order": it is that things which do not compete for
the same physical position must not be resolved in the same chain. Duricrusts are
pedogenic horizons that form WITHIN a profile and become a surface only where
erosion strips what lay above; loess and diatomite sit ON the profile. One
precedence list forces an arbitrary answer to a question that has none.

So precedence exists only within an axis, is declared in
`pedology/config/surface_classes.yaml` as an ordered list, and every region
records which rule fired AND which other rules of that axis it also satisfied.
A rule that never fires and a rule that always fires are both bugs, and neither
is visible without that record.

## There is no wet/dry axis, and there is not going to be one here

`surface_cover` carries INUNDATION, as `water`, from the solved lake extent.
It does not carry whether the ground under an unflooded class is wet or dry, and
that is a decision rather than a gap: a wet/dry cover state would be a
near-surface water state defined by this component, which is the second central
hydrology DUST-17 exists to prevent. The one authoritative top-layer liquid
water state is DUST-17's, off the land column, and a consumer that needs the
wetness of a surface reads that at its own declared depth. The class table's job
is what the surface IS MADE OF; how much water is in it at a given moment is the
column's.

The consumer this most affects is the modelled soil albedo, which is constant in
time for related reasons. `pedology/config/surface_classes.yaml` states the
decision beside the `water` rule and
`exoplasim/notes/soil-albedo-moisture.md` carries the argument.

## What can fail here

Seven checks, all with a right answer:

- **partition.** Every land region gets exactly one value on each axis and the
  class areas sum to the land area. A conservation law.
- **the strandline against the hypsometry.** Diatomite and pan-margin silcrete
  sit in the band between the solved lake surface and the spill level, summed
  here over mesh regions by `filled_km`. Hydrography reaches the same area from
  its own per-basin hypsometry curve, by a different field and a different code
  path, so a datum mismatch between the two elevation conventions or a sign error
  in the band shows up here and nowhere else.
- **the empty slot.** `surface_class == 2` is INLAND_WATER and Orogen leaves it
  deliberately empty, which is the reason this component owns standing water.
  If a future build stops leaving it empty, the water rule is no longer the only
  writer and this says so.
- **bare and soil against the albedo split.** `bare` reads the barren rock
  classes from `config/planet.yaml`, which is the same list
  `build_surface_albedo.py` divides land on. The check is stated on the two
  fallback classes rather than on their areas: no `soil` region may sit on a
  barren lithology and no `bare` region off one. Water, diatomite and pavement
  legitimately take barren cells ahead of `bare` -- a closed basin's strandline
  IS playa, and calling it diatomite is what the class is for -- so an area
  identity would have been the wrong statement, and saying so is what the first
  run of this check established.
- **gypcrete is drier than calcrete.** Watson (1983): gypsum crusts occupy the
  driest zones, calcretes the less arid parts. An inverted monthly test would
  reverse it.
- **gypcrete is enriched in carbonate-poor catchments.** The design predicts
  gypcrete rare and concentrated where the catchment already drains evaporite,
  because the alkaline path removes Ca at the first calcite. If the ion gate is
  doing nothing the enrichment comes out at 1.
- **the rule audit.** No rule fires on none of the land or almost all of it
  without being reported as such.

Comparing this classifier against another classifier would only ever tell us they
differ. `docs/src/practice/failure-modes.md` class 17.
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

from _paths import (ANALYSIS, COMPONENT_ROOT, CONFIG, PROJECT_ROOT,
                    climatology_path)

import climatology as climatology_lib  # noqa: E402  from lib/, via _paths

import builds
from gridding import climatology_cells
from orogen import INLAND_WATER, LAND, Export
from paths import rel
from provenance import require_build, staged_surface_field

from lithology_map import MEYBECK, ROCK_TO_MEYBECK
import solute_routing
from solute_routing import SPECIES

RULES = COMPONENT_ROOT / "config" / "surface_classes.yaml"

# Enumerations written into the netCDF. Order is the declared precedence for the
# rules, with the fallbacks last; the integer codes are the file's contract and
# should not be renumbered without a schema bump.
COVER = {"water": 0, "loess": 1, "diatomite": 2, "pavement": 3, "bare": 4,
         "soil": 5, "not_land": -1}
DURICRUST = {"none": 0, "gypcrete": 1, "calcrete": 2, "silcrete": 3,
             "not_land": -1}


def load_penman():
    """Hydrography's Penman implementation, imported deliberately.

    A cross-component import and it needs care, because both components ship a
    private `_paths` and a plain import would hand hydrography's module
    pedology's ANALYSIS and DATA. The swap below is explicit and is reversed in a
    `finally`.

    Reusing it rather than writing a second Penman is the point. This one is
    checked against the model's own evaporation over ocean cells, which ARE open
    water, so it has a right answer somewhere; a second implementation here would
    have none and the two could only differ.
    """
    hydro = PROJECT_ROOT / "hydrography" / "scripts"
    saved = sys.modules.pop("_paths", None)
    sys.path.insert(0, str(hydro))
    try:
        spec = importlib.util.spec_from_file_location(
            "_carve_verdict", hydro / "carve_verdict.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(hydro))
        sys.modules.pop("_paths", None)
        if saved is not None:
            sys.modules["_paths"] = saved
    return module


def monthly_climate(climatology: Path, config: dict, gravity: float):
    """Per-bin precipitation and potential evaporation, mm per Earth year.

    Vesper's year is 183 days in 12 output bins, so a bin is 15.2 days and the
    twelve-way partition of the year is the analogue of the twelve months
    Watson's criterion is stated on. An annual mean would pass cells whose wet
    season disqualifies them, which is the whole reason the criterion is monthly.

    Potential evaporation is open-water Penman, the project's own, evaluated bin
    by bin on that bin's own fields rather than on annual means. `es` is convex
    in temperature at about 6.7% per kelvin, so an annual mean understates the
    vapour pressure deficit exactly as a daily mean does, and the same Jensen
    argument that put the diurnal integration into Penman applies to the season.
    """
    penman = load_penman()
    with nc.Dataset(climatology) as ds:
        centres = np.asarray(ds["time"][:], dtype=float)
        lat = np.asarray(ds["lat"][:], dtype=float)
        scale = 1000.0 * 86400.0 * solute_routing.EARTH_YEAR_DAYS
        pr = np.asarray(ds["pr"][:], dtype=float) * scale
        ts = np.asarray(ds["ts"][:], dtype=float)
        tas = np.asarray(ds["tas"][:], dtype=float)
        ps = np.asarray(ds["ps"][:], dtype=float) * 100.0
        rss = np.asarray(ds["rss"][:], dtype=float)
        rls = np.asarray(ds["rls"][:], dtype=float)
        q_air = np.asarray(ds["hus"][:], dtype=float)[:, -1]
        wind = np.asarray(ds["spd"][:], dtype=float)[:, -1]
        diurnal = (np.asarray(ds["maxt"][:], dtype=float)
                   - np.asarray(ds["mint"][:], dtype=float))

    # The staged background land albedo, resolved through the one door rather
    # than by rebuilding the path from `model.resolution`. That path is keyed by
    # the RUNG alone while `surface_albedo` rewrites it per BUILD, so it cannot
    # say which build's field is in it, and this read used to take whichever was
    # there. `staged_surface_field` refuses a field from another build and hands
    # back the record that goes into this script's own report, so the pairing is
    # auditable after the fact. world-xgtj, CLAUDE.md rule 5.
    staged_albedo = staged_surface_field(174, config)
    land_albedo = penman.read_sra_field(
        PROJECT_ROOT / staged_albedo["path"], *ps.shape[1:])

    pet = np.empty_like(pr)
    for t in range(pr.shape[0]):
        pet[t] = penman.penman_open_water(
            ts[t], tas[t], q_air[t], wind[t], ps[t], rss[t], rls[t],
            land_albedo, gravity, diurnal_range=diurnal[t]) * scale
    return lat, pr, pet, centres, staged_albedo


def region_chemistry(export: Export, terminal: np.ndarray, brine: dict,
                     n_basins: int):
    """Per-region Ca/HCO3 and dissolved silica, from the basin or from the rock.

    A cell that drains to a closed basin inherits its basin's flux-weighted
    catchment chemistry, because that is the water that arrives and evaporates
    there. A cell that drains to the ocean has no basin, so it gets the release
    of the rock beneath it: there is no evaporative concentration to inherit and
    the local supply is what a pedogenic horizon in that profile sees.

    A cell inside a basin whose catchment generates NO runoff also falls back to
    the local rock, and that fallback is reported rather than left silent: no
    water arrives, so `brine_paths.py` reports the divide undetermined there, and
    what a pedogenic horizon sees is what is beneath it. The count is in the
    report because a silent fallback across a component boundary is how a
    plausible number replaces an error.

    Returns (ca_over_hco3, silica_umol_l, local_ca_ueq_l, has_basin_chemistry).
    """
    tables = json.loads(MEYBECK.read_text(encoding="utf-8"))["table_2c"]
    cols = {c: i for i, c in enumerate(tables["columns"])}
    release = solute_routing.region_release(export, ROCK_TO_MEYBECK,
                                            tables["rows"], cols)
    ca = release[SPECIES.index("ca_ueq_l")]
    hco3 = release[SPECIES.index("hco3_ueq_l")]
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(hco3 > 0, ca / np.maximum(hco3, 1e-30), np.nan)
    silica = release[SPECIES.index("sio2_umol_l")].copy()

    basin_ratio = np.full(n_basins, np.nan)
    basin_silica = np.full(n_basins, np.nan)
    index = {b: i for i, b in enumerate(brine["_basin_order"])}
    for row in brine["basins"]:
        i = index[row["basin_id"]]
        basin_ratio[i] = row["ca_over_hco3"]
        basin_silica[i] = row["sio2_umol_l"]

    inside = terminal >= 0
    from_basin = inside & np.isfinite(basin_ratio[np.clip(terminal, 0, None)])
    ratio = np.where(from_basin, basin_ratio[np.clip(terminal, 0, None)], ratio)
    silica = np.where(from_basin, basin_silica[np.clip(terminal, 0, None)],
                      silica)
    return ratio, silica, ca, from_basin


def audit_rules(masks: dict[str, np.ndarray], area: np.ndarray,
                land: np.ndarray, limits: dict) -> dict:
    """Fired-and-satisfied per rule, with the never/always test applied."""
    total = float(area[land].sum())
    out = {}
    for name, sat in masks.items():
        share = float(area[land & sat].sum() / total)
        out[name] = {
            "land_fraction_satisfied": share,
            "never_fires": share <= limits["min_land_fraction"],
            "always_fires": share >= limits["max_land_fraction"],
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--climatology", type=Path, default=None)
    ap.add_argument("--dust", type=Path, default=None,
                    help="aeolian dust product; defaults to the baseline")
    ap.add_argument("--dust-end", choices=("low", "central", "high"),
                    default="central",
                    help="which end of the aeolian roughness bracket decides "
                         "the loess and pavement rules. `low` is the SMOOTH bed "
                         "and therefore the DUSTY end; the names are the z0 "
                         "bracket's, not the deposition's.")
    ap.add_argument("--brine", type=Path, default=None)
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()

    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    rules = yaml.safe_load(RULES.read_text(encoding="utf-8"))
    build = str(config["source_build"])
    gravity = float(config["planet"]["gravity_m_s2"])

    climatology = args.climatology or climatology_path()
    if not climatology.is_file():
        raise SystemExit(f"{climatology} does not exist")
    require_build(climatology, "climatology", config)

    dust = args.dust or (PROJECT_ROOT / "aeolian" / "analysis"
                         / "dust_baseline.nc")
    # REFUSES rather than warns, and this is the input that has to. The
    # climatology is checked, the lake solution and the drainage are namespaced
    # per build so the wrong one is not at the path at all, and
    # `brine_paths.json` states its build. The dust field is the only input
    # living at a fixed path with no build in it, so it is the one that would be
    # silently inherited from a superseded terrain. Re-running
    # `aeolian/scripts/build_dust.py` is what makes it current; there is no
    # per-build copy to fall back to.
    require_build(dust, "dust field", config, allow_unstamped=False)

    brine_path = args.brine or (ANALYSIS / "brine_paths.json")
    if not brine_path.is_file():
        raise SystemExit(
            f"{rel(brine_path)} does not exist; run brine_paths.py first. It "
            "carries the per-basin chemical divide the duricrust ion gates and "
            "the silica gates read, and there is deliberately no fallback: "
            "guessing an ion source would decide which evaporite a basin grows.")
    brine = json.loads(brine_path.read_text(encoding="utf-8"))
    if brine.get("source_build") != build:
        raise SystemExit(
            f"{rel(brine_path)} was computed on build {brine.get('source_build')!r} "
            f"and this is {build!r}")

    export = Export(builds.mesh_export(config))
    grid_dir = builds.grid_export(config)
    hydro = builds.component_data("hydrography", config, strict=True)

    land = export.surface_class == LAND
    area = export.cell_area.astype(np.float64)
    substrate = export.substrate_class.astype(np.int32)
    codes = {c["code"]: int(c["id"])
             for c in export.manifest["lithology"]["rockClasses"]}

    # The design's claim that this component owns standing water rests on Orogen
    # leaving `surface_class == 2` empty. Check it rather than inherit it.
    inland_water_regions = int((export.surface_class == INLAND_WATER).sum())

    with nc.Dataset(hydro / "regions.nc") as ds:
        terminal = np.asarray(ds["terminal"][:]).astype(np.int64)
        filled_km = np.asarray(ds["filled_km"][:], dtype=float)
    with nc.Dataset(hydro / "basins.nc") as ds:
        n_basins = ds.dimensions["basin"].size
        spill_km = np.asarray(ds["spill_km"][:], dtype=float)
        sink_elevation_km = np.asarray(ds["sink_elevation_km"][:], dtype=float)
        basin_ids = [str(x) for x in np.asarray(ds["basin_id"][:])]
    lakes_file = hydro / "surface_water.nc"
    require_build(lakes_file, "lake solution", config)
    with nc.Dataset(lakes_file) as ds:
        lake = np.asarray(ds["lake"][:]).astype(bool)
        discharge = np.asarray(ds["discharge_m3_s"][:], dtype=float)
        basin_level_km = np.asarray(ds["basin_level_km"][:], dtype=float)
        fills_to_spill = np.asarray(ds["basin_fills_to_spill"][:]).astype(bool)
    brine["_basin_order"] = basin_ids

    # --- climate, per region --------------------------------------------
    (clim_lat, pr_bins, pet_bins, bin_centres,
     staged_albedo) = monthly_climate(climatology, config, gravity)
    row, col = climatology_cells(export, grid_dir, clim_lat)
    # Weighted: the bins are not equal in length. CLIM-13.
    precip = climatology_lib.annual_mean(pr_bins, bin_centres)[row, col]
    pet_exceeds_p = (pet_bins > pr_bins).all(axis=0)[row, col]

    with nc.Dataset(dust) as ds:
        if not np.allclose(np.asarray(ds["lat"][:], dtype=float), clim_lat):
            raise SystemExit(
                "the dust field and the climatology are on different latitude "
                "rows; they must share a grid")
        deposition = {end: np.asarray(ds[f"deposition_{end}"][:], dtype=float)
                      for end in ("low", "central", "high")}

    ca_ratio, silica, local_ca, from_basin = region_chemistry(
        export, terminal, brine, n_basins)

    # --- surface_cover rules --------------------------------------------
    cover_cfg = rules["surface_cover"]
    barren_codes = config["model"]["barren_rock_classes"]
    barren = np.zeros(substrate.shape, dtype=bool)
    for code in barren_codes:
        if code in codes:
            barren |= substrate == codes[code]

    # Strandline: inside a basin, above the solved lake, below the spill. Empty
    # for a basin with no lake and empty for one pinned at its spill, so the two
    # disqualifying cases fall out of the geometry.
    alternating = (basin_level_km > sink_elevation_km) & (~fills_to_spill)
    safe_terminal = np.clip(terminal, 0, None)
    strandline = (
        (terminal >= 0) & (~lake) & alternating[safe_terminal]
        & (filled_km > basin_level_km[safe_terminal])
        & (filled_km <= spill_km[safe_terminal]))

    dep = deposition[args.dust_end][row, col]
    loess_threshold = float(cover_cfg["loess"]["deposition_g_m2_yr"])
    traps = (~barren) & (precip >= cover_cfg["loess"]
                         ["trapping_precipitation_mm_yr"])

    cover_masks = {
        "water": lake,
        "loess": (dep >= loess_threshold) & traps,
        "diatomite": strandline & (silica >= cover_cfg["diatomite"]
                                   ["minimum_silica_umol_l"]),
        "pavement": ((~barren)
                     & (precip <= cover_cfg["pavement"]
                        ["maximum_precipitation_mm_yr"])
                     & (dep < loess_threshold)),
        "bare": barren,
    }
    cover = np.full(substrate.shape, COVER["not_land"], dtype=np.int8)
    cover[land] = COVER["soil"]
    cover_satisfied = np.zeros(substrate.shape, dtype=np.uint16)
    for bit, name in enumerate(cover_cfg["precedence"]):
        if name not in cover_masks:
            continue
        cover_satisfied |= (land & cover_masks[name]).astype(np.uint16) << bit
    for name in reversed(cover_cfg["precedence"]):
        if name in cover_masks:
            cover[land & cover_masks[name]] = COVER[name]

    # --- duricrust rules -------------------------------------------------
    duri_cfg = rules["duricrust"]
    carbonate_poor = np.nan_to_num(ca_ratio, nan=0.0) >= 1.0

    land_discharge = discharge[land & (discharge > 0)]
    drainage_cut = float(np.percentile(
        land_discharge, duri_cfg["silcrete"]["drainage_line_discharge_percentile"]
    )) if land_discharge.size else np.inf

    duri_masks = {
        "gypcrete": ((precip <= duri_cfg["gypcrete"]["maximum_precipitation_mm_yr"])
                     & pet_exceeds_p
                     & (carbonate_poor
                        if duri_cfg["gypcrete"]["require_carbonate_poor"]
                        else True)),
        "calcrete": ((precip >= duri_cfg["calcrete"]["minimum_precipitation_mm_yr"])
                     & (precip <= duri_cfg["calcrete"]["maximum_precipitation_mm_yr"])
                     & (local_ca >= duri_cfg["calcrete"]["minimum_ca_ueq_l"])),
        "silcrete": ((strandline | (discharge >= drainage_cut))
                     & (silica >= duri_cfg["silcrete"]["minimum_silica_umol_l"])),
    }
    # Which silcrete setting fired, kept apart because the two are different
    # settings with the same product and the design asks whether they can be
    # separated at all.
    silcrete_setting = {
        "pan_margin": duri_masks["silcrete"] & strandline,
        "drainage_line": duri_masks["silcrete"] & (discharge >= drainage_cut),
    }

    duricrust = np.full(substrate.shape, DURICRUST["not_land"], dtype=np.int8)
    duricrust[land] = DURICRUST["none"]
    duri_satisfied = np.zeros(substrate.shape, dtype=np.uint16)
    for bit, name in enumerate(duri_cfg["precedence"]):
        duri_satisfied |= (land & duri_masks[name]).astype(np.uint16) << bit
    for name in reversed(duri_cfg["precedence"]):
        duricrust[land & duri_masks[name]] = DURICRUST[name]

    # --- checks that can fail --------------------------------------------
    land_area = float(area[land].sum())

    def share(mask):
        return float(area[land & mask].sum() / land_area)

    cover_areas = {name: share(cover == COVER[name])
                   for name in cover_cfg["precedence"]}
    duri_areas = {name: share(duricrust == DURICRUST[name])
                  for name in ["none"] + list(duri_cfg["precedence"])}

    gyp = land & (duricrust == DURICRUST["gypcrete"])
    cal = land & (duricrust == DURICRUST["calcrete"])
    p_gyp = (float(np.average(precip[gyp], weights=area[gyp]))
             if gyp.any() else float("nan"))
    p_cal = (float(np.average(precip[cal], weights=area[cal]))
             if cal.any() else float("nan"))

    endorheic = land & (terminal >= 0)
    ca_rich_share = (float(area[endorheic & carbonate_poor].sum()
                           / max(area[endorheic].sum(), 1e-30)))
    gyp_in_ca_rich = (float(area[gyp & carbonate_poor].sum()
                            / max(area[gyp].sum(), 1e-30)) if gyp.any()
                      else float("nan"))

    # The strandline band against the hypsometry it ought to be the complement
    # of. Two independent routes to the same area: this one sums `cell_area` over
    # regions whose `filled_km` lies between the solved lake surface and the
    # spill, and hydrography's `area_at_spill_km2` comes off the per-basin
    # hypsometry curve rebuilt on the finished terrain. They share a terrain and
    # nothing else, so a datum mismatch between `filled_km` and `level_km`, or a
    # sign error in the band, shows up here and nowhere else.
    #
    # The tolerance is 10% on the total and is set from mesh granularity rather
    # than from the answer: a region is about 230 km2 and the band is a thin
    # annulus, so a region's width around its edge is the resolution of the
    # comparison. That is an order of magnitude looser than what this build
    # actually returns.
    band_area = np.bincount(terminal[strandline], weights=area[strandline],
                            minlength=n_basins)
    lake_in_basin = land & lake & (terminal >= 0)
    lake_area = np.bincount(terminal[lake_in_basin],
                            weights=area[lake_in_basin], minlength=n_basins)
    with nc.Dataset(hydro / "basins.nc") as ds:
        area_at_spill = np.asarray(ds["area_at_spill_km2"][:], dtype=float)
    expected_band = np.maximum(area_at_spill - lake_area, 0.0)
    banded = alternating & (expected_band > 0)
    band_ratio = (float(band_area[banded].sum() / expected_band[banded].sum())
                  if banded.any() and expected_band[banded].sum() > 0
                  else float("nan"))

    checks = {
        "strandline_matches_the_hypsometry": {
            "criterion": "the band between the solved lake surface and the "
                         "spill, summed over mesh regions, against "
                         "area_at_spill_km2 less the lake area, which comes off "
                         "hydrography's own hypsometry curve. Ratio within 10%, "
                         "a tolerance set from the 230 km2 region against a thin "
                         "annulus and not from the answer",
            "basins_tested": int(banded.sum()),
            "band_area_km2": float(band_area[banded].sum()),
            "expected_area_km2": float(expected_band[banded].sum()),
            "ratio": band_ratio,
            "basins_with_an_empty_band_despite_a_positive_expectation": int(
                (band_area[banded] == 0).sum()),
            "passes": bool(abs(band_ratio - 1.0) < 0.10),
        },
        "partition_cover": {
            "criterion": "class areas sum to the land area within 1e-9 relative",
            "residual": abs(sum(cover_areas.values()) - 1.0),
            "passes": abs(sum(cover_areas.values()) - 1.0) < 1e-9,
        },
        "partition_duricrust": {
            "criterion": "class areas sum to the land area within 1e-9 relative",
            "residual": abs(sum(duri_areas.values()) - 1.0),
            "passes": abs(sum(duri_areas.values()) - 1.0) < 1e-9,
        },
        "inland_water_slot_is_empty": {
            "criterion": "Orogen leaves surface_class == 2 empty, so standing "
                         "water is this component's to write",
            "regions": inland_water_regions,
            "passes": inland_water_regions == 0,
        },
        "bare_and_soil_respect_the_albedo_split": {
            "criterion": "`bare` and `soil` are the two ends of the same split "
                         "build_surface_albedo.py makes, so no `soil` region may "
                         "sit on a barren lithology and no `bare` region may sit "
                         "off one. Stated on the two fallback classes rather than "
                         "on their areas, because water, diatomite and the rest "
                         "legitimately take barren cells ahead of `bare` -- the "
                         "strandline of a closed basin IS playa, and calling it "
                         "diatomite is the point of the class",
            "soil_on_barren_regions": int((land & (cover == COVER["soil"])
                                           & barren).sum()),
            "bare_off_barren_regions": int((land & (cover == COVER["bare"])
                                            & ~barren).sum()),
            "barren_land_fraction": share(barren),
            "barren_taken_by": {
                name: share(barren & (cover == COVER[name]))
                for name in cover_cfg["precedence"]},
            "passes": bool(not ((land & (cover == COVER["soil"]) & barren).any()
                                or (land & (cover == COVER["bare"])
                                    & ~barren).any())),
        },
        "gypcrete_is_drier_than_calcrete": {
            "criterion": "Watson (1983): gypsum crusts occupy the driest zones, "
                         "calcretes the less arid parts. Area-weighted mean "
                         "annual precipitation, gypcrete strictly below calcrete",
            "gypcrete_precip_mm_yr": p_gyp,
            "calcrete_precip_mm_yr": p_cal,
            "passes": bool(p_gyp < p_cal) if gyp.any() and cal.any() else None,
        },
        "gypcrete_is_enriched_in_carbonate_poor_catchments": {
            "criterion": "the alkaline path removes Ca at the first calcite, so "
                         "gypcrete should concentrate where the catchment is "
                         "carbonate-poor. Enrichment above 1 means the ion gate "
                         "binds; 1 means it does nothing",
            "carbonate_poor_share_of_endorheic_land": ca_rich_share,
            "carbonate_poor_share_of_gypcrete": gyp_in_ca_rich,
            "enrichment": (gyp_in_ca_rich / ca_rich_share
                           if ca_rich_share > 0 else None),
            "passes": (bool(gyp_in_ca_rich > ca_rich_share)
                       if gyp.any() and ca_rich_share > 0 else None),
        },
    }

    audit = {
        "surface_cover": audit_rules(cover_masks, area, land, rules["audit"]),
        "duricrust": audit_rules(duri_masks, area, land, rules["audit"]),
    }

    # Does each ion gate bind? A gate that never excludes anything is a rule that
    # always fires wearing a different shape, and the design asks the question
    # explicitly for calcrete: Ca and HCO3 are supplied by every silicate and
    # carbonate lithology in Meybeck's table, so the gate is expected NOT to bind
    # and a lithology that failed it would have to be visible rather than
    # silently absent.
    gypcrete_climate = (
        (precip <= duri_cfg["gypcrete"]["maximum_precipitation_mm_yr"])
        & pet_exceeds_p)
    calcrete_climate = (
        (precip >= duri_cfg["calcrete"]["minimum_precipitation_mm_yr"])
        & (precip <= duri_cfg["calcrete"]["maximum_precipitation_mm_yr"]))
    silcrete_setting_any = strandline | (discharge >= drainage_cut)
    gates = {
        "gypcrete_carbonate_poor": {
            "land_fraction_before_gate": share(gypcrete_climate),
            "land_fraction_after_gate": share(duri_masks["gypcrete"]),
        },
        "calcrete_ca_source": {
            "land_fraction_before_gate": share(calcrete_climate),
            "land_fraction_after_gate": share(duri_masks["calcrete"]),
        },
        "silcrete_silica_supply": {
            "land_fraction_before_gate": share(silcrete_setting_any),
            "land_fraction_after_gate": share(duri_masks["silcrete"]),
        },
    }
    for g in gates.values():
        before = g["land_fraction_before_gate"]
        g["binds"] = bool(before > 0 and
                          g["land_fraction_after_gate"] < 0.999 * before)

    # How much of the silcrete SETTING a climate-keyed crust takes ahead of it.
    # Recoverable from the bitmask, and recorded because the ordering is a
    # declared judgment rather than a fact and a reader should see its cost.
    silcrete_under = {
        name: share(duri_masks["silcrete"] & (duricrust == DURICRUST[name]))
        for name in ["gypcrete", "calcrete"]}

    # Loess and pavement across the aeolian roughness bracket, because that
    # bracket is the widest thing in the dust product and quoting the central
    # value alone would hide a factor of 40 in emission.
    bracket = {}
    for end, dep_end in deposition.items():
        d = dep_end[row, col]
        bracket[end] = {
            "loess_land_fraction": share((d >= loess_threshold) & traps),
            "pavement_land_fraction": share(
                (~barren)
                & (precip <= cover_cfg["pavement"]["maximum_precipitation_mm_yr"])
                & (d < loess_threshold)),
            "land_mean_deposition_g_m2_yr": float(
                np.average(d[land], weights=area[land])),
        }

    # --- write ------------------------------------------------------------
    out = args.output or (COMPONENT_ROOT / "analysis" / "surface_classes.nc")
    out.parent.mkdir(parents=True, exist_ok=True)
    with nc.Dataset(out, "w", format="NETCDF4") as ds:
        ds.createDimension("region", substrate.size)
        v = ds.createVariable("surface_cover", "i1", ("region",), zlib=True)
        v[:] = cover
        v.long_name = "what wind and light see"
        v.flag_values = np.array([COVER[k] for k in COVER], dtype=np.int8)
        v.flag_meanings = " ".join(COVER)
        v = ds.createVariable("duricrust", "i1", ("region",), zlib=True)
        v[:] = duricrust
        v.long_name = "what is cementing at or below the surface"
        v.flag_values = np.array([DURICRUST[k] for k in DURICRUST], dtype=np.int8)
        v.flag_meanings = " ".join(DURICRUST)
        v = ds.createVariable("surface_cover_satisfied", "u2", ("region",),
                              zlib=True)
        v[:] = cover_satisfied
        v.long_name = ("bitmask, bit i set where rule i of the surface_cover "
                       "precedence was satisfied whether or not it fired")
        v.bit_meanings = " ".join(cover_cfg["precedence"])
        v = ds.createVariable("duricrust_satisfied", "u2", ("region",),
                              zlib=True)
        v[:] = duri_satisfied
        v.long_name = ("bitmask, bit i set where rule i of the duricrust "
                       "precedence was satisfied whether or not it fired")
        v.bit_meanings = " ".join(duri_cfg["precedence"])
        setting = np.zeros(substrate.shape, dtype=np.uint16)
        for bit, name in enumerate(duri_cfg["silcrete"]["settings"]):
            setting |= (land & silcrete_setting[name]).astype(np.uint16) << bit
        v = ds.createVariable("silcrete_setting", "u2", ("region",), zlib=True)
        v[:] = setting
        v.long_name = ("bitmask of which non-pedogenic silcrete setting is "
                       "present; the types share one duricrust value because it "
                       "is not clear this resolution separates their products, "
                       "and this keeps the distinction that does exist")
        v.bit_meanings = " ".join(duri_cfg["silcrete"]["settings"])
        v.groundwater_note = ("Ullyott and Nash's groundwater type is NOT placed: "
                              "it sits at or near a water table and this project "
                              "models none. The slot is empty rather than filled "
                              "with a surface-elevation proxy.")
        ds.setncattr("vesper_source_build", build)
        ds.setncattr("vesper_terrain_hash", export.terrain_hash)
        ds.setncattr("climatology", rel(climatology))
        ds.setncattr("dust_field", rel(dust))
        ds.setncattr("dust_bracket_end", args.dust_end)
        ds.setncattr("note", "Derived surface classes on two independent axes. "
                             "Computed from a PRE-CARVE climatology and lake "
                             "solution; regenerate after the next terrain "
                             "iteration.")

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": "pedology/scripts/build_surface_classes.py",
        "source_build": build,
        "terrain_hash": export.terrain_hash,
        "climatology": rel(climatology),
        "climatology_sha256": hashlib.sha256(climatology.read_bytes()).hexdigest(),
        # Which build's staged background albedo the Penman evaporation was
        # taken over. The path is keyed by the rung alone, so recording it by
        # name would say nothing; the record carries the terrain hash and the
        # file's own sha256. world-xgtj.
        "staged_background_albedo": staged_albedo,
        "dust_field": rel(dust),
        "dust_bracket_end": args.dust_end,
        "brine_paths": rel(brine_path),
        "rules_sha256": hashlib.sha256(RULES.read_bytes()).hexdigest(),
        "git_commit": subprocess.run(["git", "rev-parse", "HEAD"],
                                     capture_output=True, text=True,
                                     cwd=PROJECT_ROOT).stdout.strip() or None,
        "land_area_km2": land_area,
        "surface_cover_land_fractions": cover_areas,
        "duricrust_land_fractions": duri_areas,
        "silcrete_settings_land_fractions": {
            k: share(v) for k, v in silcrete_setting.items()},
        "drainage_line_discharge_cut_m3_s": drainage_cut,
        "chemistry_source": {
            "note": "A region inside a closed basin inherits its basin's "
                    "flux-weighted catchment chemistry, because that is the "
                    "water that arrives and evaporates there. A region draining "
                    "to the ocean, and a region in a basin whose catchment "
                    "generates no runoff at all, take the release of the rock "
                    "beneath them instead: there is no evaporative "
                    "concentration to inherit. The fallback is counted here "
                    "rather than left silent.",
            "land_fraction_from_basin": share(from_basin),
            "land_fraction_from_local_rock": share(~from_basin),
            "land_fraction_in_a_basin_with_no_runoff": share(
                (terminal >= 0) & ~from_basin),
        },
        "checks": checks,
        "rule_audit": audit,
        "ion_and_supply_gates": gates,
        "silcrete_setting_taken_by_a_climate_keyed_crust": silcrete_under,
        "aeolian_roughness_bracket": bracket,
        "disposable": ("Every fraction here is computed from a pre-carve "
                       "climatology, lake solution and dust field. The carve "
                       "replaces all three, so these are regenerated rather "
                       "than migrated. What survives is the classifier."),
    }
    (ANALYSIS / "surface_classes_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"{'surface_cover':>16}  {'% land':>7}")
    for name, frac in cover_areas.items():
        print(f"{name:>16}  {100*frac:7.3f}")
    print(f"\n{'duricrust':>16}  {'% land':>7}")
    for name, frac in duri_areas.items():
        print(f"{name:>16}  {100*frac:7.3f}")
    print("\nchecks")
    for name, c in checks.items():
        state = {True: "PASS", False: "FAIL", None: "n/a"}[c["passes"]]
        print(f"  {name:48} {state}")
    print("\nloess and pavement across the aeolian roughness bracket")
    for end, b in bracket.items():
        print(f"  z0 {end:>7}: deposition {b['land_mean_deposition_g_m2_yr']:9.3f} "
              f"g/m2/yr  loess {100*b['loess_land_fraction']:6.3f}%  "
              f"pavement {100*b['pavement_land_fraction']:6.3f}%")
    print(f"\nwrote {rel(out)}")
    print(f"wrote {rel(ANALYSIS / 'surface_classes_report.json')}")


if __name__ == "__main__":
    main()
