#!/usr/bin/env python3
"""Absolute phosphorus release from rock, the root-zone stock, and the rank map.

    python pedology/scripts/phosphorus_budget.py

Writes `pedology/analysis/phosphorus_budget.json`.

This is a worldbuilding project: everything below is the simulated phosphorus
supply of Vesper, an invented super-Earth, computed from its lithology and its
modelled climate.

## Two arms, and the units are the difference

**The absolute arm** is a mass flux. Hartmann et al. (2014) define the chemical
weathering rate as the fluvial export of Ca, Mg, Na, K and SiO2 and give the
phosphorus release as a fixed percentage of it per lithological class, so

    F_P = release_percent/100 * F_(Ca + Mg + Na + K + SiO2)

is the whole law, and the only thing it needs that this project did not already
have is an absolute major-element flux. Meybeck (1987) Table 2C supplies the
per-lithology concentrations and runoff carries them, which is the same model
form `weathering_fluxes.py` already uses for silica and bicarbonate and the same
extracted table. Nothing is fitted, nothing is normalised, and every factor has
a unit. The root-zone stock sits beside it: parent P content times bulk density
times the regolith depth cut at the land column property contract's rootable
base, in kgP/m2. It is a TOTAL and an upper bound, not a labile pool.

**The rank arm** is what this script used to be and is kept, deliberately, as a
rank. It scores where the outbound leg may matter and where a closed basin
concentrates what the catchment delivered, per mesh region and with the drainage
terminal attached. It carries no mass and no time, and converting it into one is
the thing the absolute arm exists to make unnecessary.

The two live on different supports and are labelled with them: the absolute arm
on the climatology's grid, because runoff is a climate field, and the rank arm on
the native mesh, because drainage terminals are a mesh product. Neither is
matched to the other by longitude.

## What the absolute arm does NOT carry

Its `not_carried` block is the register, and it is a register rather than a set
of defaults because a zero and an absence are different claims. The largest
entries are a regolith production RATE, which this project has no time axis for,
a soil-shielding function, which Hartmann applies and no source here supplies,
and the sorption and retention that stand between release and a root.

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
from _paths import ANALYSIS, CONFIG, PROJECT_ROOT, climatology_path  # noqa: E402
from paths import rel  # noqa: E402

sys.path.insert(0, str(PROJECT_ROOT / "lib"))
import builds  # noqa: E402
from orogen import Export, LAND  # noqa: E402
from surface_classes import cover_mask  # noqa: E402

import yaml  # noqa: E402

PEDOGENESIS = HERE.parent / "config" / "pedogenesis.yaml"
MEYBECK = PROJECT_ROOT / "pedology" / "data" / "reference" / "meybeck1987_tables.json"
LAND_COLUMN = HERE.parent / "config" / "land_column_properties.yaml"
OUT = ANALYSIS / "phosphorus_budget.json"

# --------------------------------------------------------------------------
# The absolute arm's constants, and the two criteria it is judged against.
# Both criteria are fixed here, before the arm runs, because a bound chosen
# after the result it judges is not a bound.
# --------------------------------------------------------------------------

# Hartmann's weathering rate is the fluvial export of Ca, Mg, Na, K and SiO2.
# Meybeck's Table 2C carries the cations in ueq/l and silica in umol/l, which is
# why they are converted here per species rather than summed through the table's
# own `cation_sum_ueq_l`: an equivalent is not a mole and not a gram.
CATIONS = {                      # Meybeck column -> (g/mol, charge)
    "ca_ueq_l": (40.078, 2),
    "mg_ueq_l": (24.305, 2),
    "na_ueq_l": (22.990, 1),
    "k_ueq_l": (39.098, 1),
}
SIO2_G_PER_MOL = 60.084

# Hartmann et al. (2014) section 2.1 excludes evaporite dissolution beyond
# carbonate from the weathering rate. Meybeck's two evaporite rows are almost
# entirely halite and gypsum back in solution, so their cation load is salt
# dissolving rather than rock weathering. For those two classes the cation flux
# is cut to the share balanced by BICARBONATE, which is the part hydrolysis
# produced, and the cut is applied to every cation in proportion to its own
# equivalents so that no species split is invented.
EVAPORITE_CLASSES = {"halite_evaporite", "gypsum_evaporite"}

# CRITERION 1, declared before the arm ran. Hartmann and Moosdorf (2011) model
# P release on the Japanese Archipelago -- wet, warm and volcanic, the high end
# of Earth's range -- and reach 390 kgP/km2/yr. This world's land is drier than
# Earth's. A land mean above that figure would be claiming Vesper out-releases
# the wettest volcanic archipelago on Earth, and is reported as a failure rather
# than as a result.
HARTMANN_JAPAN_MAX_KG_KM2_YR = 390.0

# CRITERION 2, also declared first. `phosphorus_ppm` is the rock's P content and
# `phosphorus_release_relative` is that content as a PERCENTAGE of the rock's
# Ca+Mg+Na+K+SiO2 content, so the two rows over-determine the implied major
# element content of each class: content_ppm * 1e-4 / release_percent. Every
# class the map carries must land inside this range, or the two rows are not the
# same paper's quantities and the arm is not entitled to multiply one by a flux
# of the other. This is the check that the pair which has been confused twice in
# this file's history is still the pair it claims to be.
IMPLIED_MAJOR_FRACTION_RANGE = (0.40, 0.95)

# Earth's own answer to the same question, per lithological class, from the same
# paper the law comes from: Hartmann et al. (2014) Table 3, averaged P-release
# rate in kgP/km2/yr, global sum 1144e6 kgP/yr over 134e6 km2 for a global mean
# of 8.53. A COMPARATOR AND NOT A CRITERION, and it was added after the arm
# first ran, which is exactly why it is not one: the two criteria above were
# fixed before. It is also not an identity to be reproduced. Hartmann's rates
# carry Earth's climate and his soil-shielding term, and this arm has neither,
# so a class that differs is informative and a class that differs by orders is
# a reason to look.
HARTMANN_TABLE3_KG_KM2_YR = {
    "granite": 6.63, "granodiorite": 10.92,          # PA, PI
    "morb": 29.77, "oib": 29.77, "flood_basalt": 29.77, "arc_basalt": 29.77,
    "arc_andesite": 37.54, "rift_bimodal": 37.54,    # VB, VI
    "gneiss": 4.41, "schist": 4.41, "quartzite": 4.41,        # MT
    "melange": 6.56,                                          # SM
    "shelf_clastic": 4.48, "foreland_clastic": 4.48,
    "continental_clastic": 4.48,                              # SS
    "carbonate": 14.46,                                       # SC
    "playa_clastic": 6.28,                                    # SU
    "evaporite": 1.49,                                        # EV
}
HARTMANN_TABLE3_GLOBAL_MEAN = 8.53


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_commit() -> str | None:
    import subprocess
    out = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                         text=True, cwd=PROJECT_ROOT).stdout.strip()
    return out or None


# --------------------------------------------------------------------------
# The absolute arm: rock to solution, and the root-zone stock
# --------------------------------------------------------------------------

def major_element_kg_per_litre() -> dict[str, float]:
    """Meybeck Table 2C as kg of Ca+Mg+Na+K+SiO2 per litre, per rock class.

    This is the DENOMINATOR of Hartmann's phosphorus law. Hartmann et al.
    (2014) define the weathering rate as the fluvial export of the major
    cations plus silica and give the P release as a fixed percentage of it, so
    an absolute P flux needs an absolute flux of exactly this quantity and of
    nothing else. Meybeck's Table 2C supplies it per lithology, which is why the
    two papers compose: the classes are the same classes, and this project
    already maps its own lithology onto both.

    Units are the trap the table's own caution line names. Silica is umol/l and
    every other species is ueq/l, so an equivalent is divided by its charge
    before it is multiplied by a molar mass. Summing the table's
    `cation_sum_ueq_l` instead would weigh a calcium ion the same as a sodium
    one.
    """
    table = json.loads(MEYBECK.read_text(encoding="utf-8"))["table_2c"]
    columns = table["columns"]
    index = {name: columns.index(name) for name in columns}
    out: dict[str, float] = {}
    for klass, row in table["rows"].items():
        # umol/l of each cation, and of silica as SiO2.
        moles = {name: row[index[name]] / charge
                 for name, (_, charge) in CATIONS.items()}
        if klass in EVAPORITE_CLASSES:
            # Salt going back into solution is not rock being weathered. Keep
            # the share of the cation load that bicarbonate balances, which is
            # the hydrolysis part, and cut every cation by the same factor.
            equivalents = sum(row[index[name]] for name in CATIONS)
            keep = (min(1.0, row[index["hco3_ueq_l"]] / equivalents)
                    if equivalents > 0 else 0.0)
            moles = {name: value * keep for name, value in moles.items()}
        grams_per_litre = sum(moles[name] * mass * 1e-6
                              for name, (mass, _) in CATIONS.items())
        grams_per_litre += row[index["sio2_umol_l"]] * SIO2_G_PER_MOL * 1e-6
        out[klass] = grams_per_litre * 1e-3          # g/l -> kg/l
    return out


def implied_major_fraction(p_ppm: dict, p_rel: dict, classes) -> dict:
    """CRITERION 2. What the two configured rows say each rock is made of.

    `phosphorus_ppm` is a content and `phosphorus_release_relative` is that
    content as a percentage of the rock's Ca+Mg+Na+K+SiO2 content, so their
    ratio is the major-element content itself and is a property of the rock
    rather than of this world. It is checked rather than reported because the
    two rows have been confused for one another twice in this file's history,
    in opposite directions, and the absolute arm multiplies one of them by a
    flux of the other. A class outside the range is a refusal.
    """
    lo, hi = IMPLIED_MAJOR_FRACTION_RANGE
    implied, outside = {}, []
    for code in sorted(classes):
        release = float(p_rel.get(code, 0.0))
        if release <= 0.0:
            continue
        value = float(p_ppm.get(code, 0.0)) * 1e-4 / release
        implied[code] = round(value, 4)
        if not lo <= value <= hi:
            outside.append(f"{code} implies {value:.3f}")
    return {"range": list(IMPLIED_MAJOR_FRACTION_RANGE),
            "by_class": implied, "outside_range": outside,
            "what_it_is": "content_ppm * 1e-4 / release_percent, which is the "
                          "mass fraction of the rock that is Ca, Mg, Na, K and "
                          "SiO2 if the two configured rows are the same paper's "
                          "quantities"}


def root_zone_stock(soil_map: Path, content_ppm: np.ndarray, land: np.ndarray,
                    lon: np.ndarray, lat: np.ndarray, decimals: int,
                    rootable_base_m: float) -> dict:
    """Total parent-derived phosphorus in the rooted column, kgP/m2.

    Depth and bulk density come from the soil map, the rooted base from the land
    column property contract, and the concentration from the lithology. It is
    the stock a root zone could draw on if every gram of the parent's phosphorus
    were still there and still reachable, which is an UPPER BOUND in both
    directions that matter: weathering has been exporting phosphorus out of the
    profile for as long as the profile has existed, and only a small share of
    what remains is ever labile. The labile state is `anut-10`'s and this is
    deliberately not it.

    The soil map is keyed by rounded coordinates, the way every other consumer
    of it is keyed, and a land cell with no row is a refusal rather than a hole:
    the two files would then describe different land.
    """
    rows = np.genfromtxt(soil_map, names=True)
    key = {(round(float(a), decimals), round(float(b), decimals)): i
           for i, (a, b) in enumerate(zip(rows["Lon"], rows["Lat"]))}
    signed = np.where(lon > 180.0, lon - 360.0, lon)

    js, is_ = np.nonzero(land)
    depth = np.empty(js.size)
    density = np.empty(js.size)
    for n, (j, i) in enumerate(zip(js, is_)):
        found = key.get((round(float(signed[i]), decimals),
                         round(float(lat[j]), decimals)))
        if found is None:
            raise SystemExit(
                f"{soil_map} has no row for the land cell at "
                f"{signed[i]:.4f}, {lat[j]:.4f}. The soil map and the "
                "climatology describe different land; rebuild the soil map "
                "against this climatology.")
        depth[n] = rows["depth"][found]
        density[n] = rows["bulkdensity"][found]

    rooted_m = np.minimum(depth, rootable_base_m)
    stock = content_ppm[land] * 1e-6 * density * rooted_m      # kgP/m2
    return {
        "definition": "parent P content times bulk density times the rooted "
                      "thickness, which is the regolith depth cut at the land "
                      "column property contract's rootable base",
        "rootable_base_m": rootable_base_m,
        "is_an_upper_bound": (
            "TOTAL phosphorus, not labile and not available. It assumes the "
            "regolith still holds its parent's P content, and weathering has "
            "been exporting P out of the profile throughout its existence. The "
            "labile pool a vegetation model consumes is a small share of this "
            "and is anut-10's to define."),
        "kg_p_per_m2": {
            "p10": float(np.percentile(stock, 10)),
            "median": float(np.median(stock)),
            "p90": float(np.percentile(stock, 90)),
            "mean": float(stock.mean()),
        },
        "rooted_thickness_m_median": float(np.median(rooted_m)),
        "cells_where_regolith_is_shallower_than_the_rootable_base": int(
            np.count_nonzero(depth < rootable_base_m)),
        "land_cells": int(js.size),
    }


def absolute_arm(cfg: dict, ped: dict, climatology: Path,
                 soil_map: Path) -> dict:
    """Rock to solution in kgP per year, and the root-zone stock beside it.

    THE LAW IS HARTMANN'S OWN, and it is one multiplication:

        F_P = b_i * F_(Ca + Mg + Na + K + SiO2)

    with `b_i` the per-lithology release percentage `pedogenesis.yaml` already
    carries from Table A1-2 and the flux from Meybeck's Table 2C concentrations
    times this world's runoff, which is Meybeck's own model form and is what
    `weathering_fluxes.py` already uses for silica and bicarbonate. Nothing here
    is fitted and nothing is normalised: every factor is a number with a unit.

    Two supports meet in this script and they are kept apart by name. This arm
    is on the CLIMATOLOGY's grid, because runoff is what carries phosphorus into
    solution and runoff is a climate field. The relative arm below it is on the
    native MESH, because drainage terminals are a mesh product. Neither is
    matched to the other by longitude; the lithology reaches this grid through
    `gridding.land_fraction_of_class`, which bins by index.
    """
    from provenance import require_build
    from climatology import annual_mean
    from gridding import land_fraction_of_class
    from build_soil import (COORD_DECIMALS, EARTH_YEAR_DAYS, KELVIN,
                            lithology_fractions)
    from weathering_fluxes import grid_cell_area_km2
    from lithology_map import ROCK_TO_MEYBECK

    require_build(climatology, "climatology", cfg)
    with Dataset(climatology) as data:
        lat = np.asarray(data["lat"][:], dtype=float)
        lon = np.asarray(data["lon"][:], dtype=float)
        land = np.asarray(data["lsm"][0], dtype=float) > 0.5
        centres = np.asarray(data["time"][:], dtype=float)
        temperature = annual_mean(np.asarray(data["tas"][:], dtype=float),
                                  centres) - KELVIN
        scale = 1000.0 * 86400.0 * EARTH_YEAR_DAYS
        evaporation = -annual_mean(np.asarray(data["evap"][:], dtype=float),
                                   centres)
        runoff = np.maximum(
            (annual_mean(np.asarray(data["pr"][:], dtype=float), centres)
             - evaporation) * scale, 0.0)

    # The configured rung's grid, for the reason weathering_fluxes.py states:
    # this step reads a climatology config has already pinned to that rung.
    grid_dir = builds.grid_export(cfg)
    fractions, mesh = lithology_fractions(cfg, grid_dir)
    area_km2 = grid_cell_area_km2(lat, lon,
                                  float(mesh.manifest["planet"]["radiusKm"]))
    # Litres of runoff per cell per Earth year: mm -> m -> m3 -> l over km2.
    litres = runoff * 1e-3 * area_km2 * 1e6 * 1e3

    major_c = major_element_kg_per_litre()
    p_ppm = ped["phosphorus_ppm"]
    p_rel = ped["phosphorus_release_relative"]

    release_kg = np.zeros_like(runoff)        # kgP per cell per Earth year
    major_kg = np.zeros_like(runoff)          # kg Si+cations per cell per year
    content = np.zeros_like(runoff)           # ppm P, land-fraction weighted
    per_class: dict[str, dict] = {}
    unmapped = []
    for code, frac in fractions.items():
        content += frac * float(p_ppm.get(code, 0.0))
        klass = ROCK_TO_MEYBECK.get(code)
        if klass is None:
            unmapped.append(code)
            continue
        major = major_c[klass] * frac * litres
        p = 0.01 * float(p_rel.get(code, 0.0)) * major
        major_kg += major
        release_kg += p
        class_area = float((frac * area_km2)[land].sum())
        per_class[code] = {
            "meybeck_class": klass,
            "release_percent_of_major": float(p_rel.get(code, 0.0)),
            "content_ppm": float(p_ppm.get(code, 0.0)),
            "kg_p_per_year": float(p[land].sum()),
            "land_area_km2": class_area,
            "yield_kg_per_km2_per_year": (float(p[land].sum()) / class_area
                                          if class_area > 0 else None),
            "earth_yield_kg_per_km2_per_year": HARTMANN_TABLE3_KG_KM2_YR.get(code),
        }

    # ENDORHEIC MEANS THE DRAINAGE FATE, NOT THE BASIN FLOOR, and the two are
    # not close: the floors Orogen flags as `is_endorheic` are the depressions
    # themselves, while what matters for where released phosphorus ends up is
    # every region whose drainage chain fails to reach the ocean. On this build
    # the second is more than twice the first. The rank arm below splits by the
    # same terminal, so one word means one thing in this artifact.
    hydro = builds.component_data("hydrography", cfg, strict=True)
    with Dataset(hydro / "regions.nc") as ds:
        terminal_region = np.asarray(ds["terminal"][:])
    endorheic = land_fraction_of_class(mesh, grid_dir, terminal_region >= 0)
    land_area = float(area_km2[land].sum())
    total = float(release_kg[land].sum())
    yield_kg_km2 = np.zeros_like(release_kg)
    np.divide(release_kg, area_km2, out=yield_kg_km2, where=area_km2 > 0)
    yields = yield_kg_km2[land]

    decl = yaml.safe_load(LAND_COLUMN.read_text(encoding="utf-8"))
    stock = root_zone_stock(soil_map, content, land, lon, lat, COORD_DECIMALS,
                            float(decl["geometry"]["rootable_base_m"]))

    land_mean_yield = total / land_area
    return {
        "what_it_is": "the absolute rock-to-solution phosphorus release and the "
                      "root-zone total-P stock. Every number carries a unit; "
                      "nothing here is normalised to a land mean.",
        "law": "F_P = release_percent/100 * F_(Ca+Mg+Na+K+SiO2), Hartmann et al. "
               "(2014) section 2.1 and Table A1-2, with the major-element flux "
               "as Meybeck (1987) Table 2C concentration times runoff",
        "support": "the climatology's grid; the lithology reaches it through "
                   "gridding.land_fraction_of_class, by index and never by "
                   "longitude",
        "climatology": rel(climatology),
        "climatology_sha256": sha256(climatology),
        "soil_map": rel(soil_map),
        "soil_map_sha256": sha256(soil_map),
        "meybeck_sha256": sha256(MEYBECK),
        "evaporite_rule": (
            "Hartmann excludes evaporite dissolution beyond carbonate, so the "
            "two evaporite classes keep only the share of their cation load "
            "that bicarbonate balances. Every cation is cut by the same factor, "
            "so no species split is invented."),
        "release": {
            "kg_p_per_year": total,
            "land_yield_kg_per_km2_per_year": land_mean_yield,
            "yield_percentiles_kg_per_km2_per_year": {
                "p10": float(np.percentile(yields, 10)),
                "median": float(np.median(yields)),
                "p90": float(np.percentile(yields, 90)),
                "max": float(yields.max()),
            },
            "major_element_flux_kg_per_year": float(major_kg[land].sum()),
            "endorheic_kg_p_per_year": float((release_kg * endorheic)[land].sum()),
            "endorheic_means": "the drainage chain does not reach the ocean, "
                               "from hydrography's regions.nc terminal, and NOT "
                               "the closed-basin floor area Orogen flags",
            "land_area_km2": land_area,
            "mean_temperature_c": float(temperature[land].mean()),
            "mean_runoff_mm_per_earth_year": float(runoff[land].mean()),
        },
        "by_rock_class": dict(sorted(per_class.items(),
                                     key=lambda kv: -kv[1]["kg_p_per_year"])),
        "rock_classes_with_no_meybeck_class": sorted(unmapped),
        "root_zone_stock": stock,
        "criteria": {
            "earth_upper_bound": {
                "declared_before_running": True,
                "bound_kg_per_km2_per_year": HARTMANN_JAPAN_MAX_KG_KM2_YR,
                "source": "Hartmann and Moosdorf (2011), modelled P release on "
                          "the Japanese Archipelago, the wet volcanic high end "
                          "of Earth's range",
                "land_yield_kg_per_km2_per_year": land_mean_yield,
                "passes": bool(land_mean_yield < HARTMANN_JAPAN_MAX_KG_KM2_YR),
                "land_cells_above_the_bound": int(
                    np.count_nonzero(yields > HARTMANN_JAPAN_MAX_KG_KM2_YR)),
            },
            "implied_major_element_fraction": implied_major_fraction(
                p_ppm, p_rel, per_class.keys()),
        },
        "earth_comparator": {
            "is_not_a_criterion": (
                "added after the arm first ran, which is why it is reported and "
                "not tested. The two criteria above were fixed before it ran. "
                "Hartmann's rates also carry Earth's climate and his soil "
                "shielding term, and this arm has neither, so agreement is not "
                "the thing being looked for."),
            "source": "Hartmann et al. (2014) Table 3",
            "earth_land_yield_kg_per_km2_per_year": HARTMANN_TABLE3_GLOBAL_MEAN,
            "earth_total_kg_p_per_year": 1144.0e6,
            "this_world_over_earth": land_mean_yield / HARTMANN_TABLE3_GLOBAL_MEAN,
            "runoff_is_the_expected_reason": (
                "this world's land is drier than Earth's, and the flux is a "
                "concentration times runoff, so a yield below Earth's is what "
                "the law predicts before any lithology contrast is considered."),
        },
        "uncertainty": {
            "note": "what the absolute numbers above are worth. Where an end "
                    "exists it is given; where none does, the entry says so and "
                    "names what would supply it, because a bracket invented to "
                    "look like one is worse than a declared absence.",
            "meybeck_is_one_analysis_per_class": (
                "Table 2C is titled 'Representative analyses used in the "
                "Temperate Stream Model' and the extracted table carries one "
                "value per species per class with no scatter beside it. So the "
                "concentration end of this flux has no published spread in any "
                "artifact here. NOT BRACKETED."),
            "hartmann_release_is_a_class_mean": (
                "the release row predicts a class mean and within-class spreads "
                "overlap, which pedogenesis.yaml already records. A cell is "
                "assigned its class's mean, so the map's contrast is real and a "
                "single cell's value is not. NOT BRACKETED."),
            "no_temperature_dependence": (
                "the flux is a concentration times runoff and Meybeck's "
                "concentrations are temperate-stream class means, so "
                "temperature enters this arm nowhere. It is reported above "
                "because it is the field the law WOULD use. Dessert et al. "
                "(2003) is the temperature-dependent alternative for the "
                "volcanic part and weathering_fluxes.py already prices the "
                "disagreement between the two on CO2; the same disagreement "
                "applies here and is not priced in phosphorus."),
            "calibrated_on_rivers_that_reach_the_ocean": {
                "what": "Meybeck's model is fluvial export from exorheic "
                        "continents, and the extracted table's own caution line "
                        "says so. On land that drains to a closed basin the "
                        "solute dissolves and then stays, so the concentration "
                        "describes what enters solution rather than what "
                        "leaves. That is the right quantity for a root zone and "
                        "the wrong one for an export, and this artifact is "
                        "reporting the first.",
                "endorheic_share_of_land": float(
                    (endorheic * area_km2)[land].sum() / land_area),
            },
        },
        "not_carried": {
            "note": "each of these is a term Hartmann's architecture or the "
                    "root-zone question needs and this pipeline does not "
                    "supply. Registered rather than defaulted: a zero and an "
                    "absence are different claims.",
            "fresh_mineral_renewal_rate": (
                "the supply of unweathered mineral into the weathering zone, "
                "which needs a regolith production RATE. build_soil.py's depth "
                "law is a steady state with no rate in it and Orogen has no "
                "time axis, so there is no length per time to be had here. "
                "docs/src/reference/no-time-axis.md"),
            "soil_residence_time": (
                "same reason. The stock and the flux are both computed, so a "
                "residence follows from dividing them, but the stock is a total "
                "and the flux is a release, so their quotient is not a "
                "residence of the same pool and is deliberately not reported."),
            "soil_shielding": (
                "Hartmann applies a soil-thickness term that shields fresh rock "
                "from water. Regolith depth exists here and the shielding "
                "FUNCTION does not: no read source in this repository gives one "
                "as a response to a soil thickness, so applying a shape would "
                "be choosing it."),
            "sorption_and_retention": (
                "what happens to released P before a root reaches it. The soil "
                "map carries the andic fixation share only, which is the "
                "excess over ordinary retention on andic ground and not a "
                "retention model. sdec-2 and sdec-4."),
            "aeolian_return": (
                "deposited mass and the soluble share of it. anut-3 owns the "
                "source-, size- and surface-resolved deposition; the relative "
                "arm below carries the source-composition hypothesis and "
                "cannot settle a sign at a destination."),
            "labile_pool": "anut-10, which owns the pedogenic initial state.",
        },
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=OUT)
    ap.add_argument("--surface-classes", type=Path, default=None,
                    help="pedology/analysis/surface_classes.nc. Given, the "
                         "budget reports the two derived classes that carry "
                         "phosphorus with opposite signs: `diatomite` as the "
                         "deflatable source and `pavement` as a sink. SURF-6")
    ap.add_argument("--climatology", type=Path, default=None,
                    help="the climatology the absolute arm takes runoff and "
                         "temperature from; defaults to the configured "
                         "baseline")
    ap.add_argument("--soil-map", type=Path, default=None,
                    help="the soil map the root-zone stock takes depth and "
                         "bulk density from; defaults to the configured "
                         "build's")
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
    climatology = (args.climatology or climatology_path()).resolve()
    if not climatology.is_file():
        raise SystemExit(f"{climatology} does not exist")
    soil_map = (args.soil_map or builds.soilmap(cfg)).resolve()
    if not soil_map.is_file():
        raise SystemExit(f"{soil_map} does not exist. Run "
                         "pedology/scripts/build_soil.py.")
    absolute = absolute_arm(cfg, ped, climatology, soil_map)

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
        "note": "Two arms. `absolute` is the rock-to-solution phosphorus release "
                "in kgP per year and the root-zone total-P stock in kgP/m2, on "
                "the climatology's grid. Everything under `relative_rank` is a "
                "RANK diagnostic on the native mesh, carries no mass and no "
                "time, and must not be converted into one. "
                "Generated by pedology/scripts/phosphorus_budget.py; do not edit.",
        "generated": datetime.now(timezone.utc).isoformat(),
        "build": cfg.get("source_build"),
        "terrain_hash": ex.terrain_hash,
        "inputs": {"pedogenesis.yaml": sha256(PEDOGENESIS),
                   "planet.yaml": sha256(CONFIG),
                   "land_column_properties.yaml": sha256(LAND_COLUMN)},
        "git_commit": git_commit(),
        "absolute": absolute,
        "relative_rank_note": (
            "the arm below ranks where the outbound leg may matter and where a "
            "closed basin concentrates it. It is a rank: `mean_release_flux_"
            "relative` is a content times a percentage and is not a flux. The "
            "absolute release above is what carries kilograms."),
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

    a = result["absolute"]
    rel_ = a["release"]
    print(f"ABSOLUTE, on {a['support'].split(';')[0]}")
    print(f"  P release             {rel_['kg_p_per_year']:.3e} kgP/yr, "
          f"{rel_['land_yield_kg_per_km2_per_year']:.2f} kgP/km2/yr at the land "
          f"mean")
    q = a["release"]["yield_percentiles_kg_per_km2_per_year"]
    print(f"    per cell            median {q['median']:.2f}, p90 {q['p90']:.2f}, "
          f"max {q['max']:.2f} kgP/km2/yr")
    print(f"    endorheic share     "
          f"{rel_['endorheic_kg_p_per_year'] / max(rel_['kg_p_per_year'], 1e-30):.1%}")
    stock = a["root_zone_stock"]
    print(f"  root-zone total P     {stock['kg_p_per_m2']['median']:.3f} kgP/m2 "
          f"median, {stock['kg_p_per_m2']['p10']:.3f} to "
          f"{stock['kg_p_per_m2']['p90']:.3f} p10-p90 (an upper bound; not labile)")
    crit = a["criteria"]["earth_upper_bound"]
    print(f"  against Earth's high end  {crit['bound_kg_per_km2_per_year']} "
          f"kgP/km2/yr: "
          + ("PASSES" if crit["passes"] else "FAILS")
          + f", {crit['land_cells_above_the_bound']} land cells above it")
    outside = a["criteria"]["implied_major_element_fraction"]["outside_range"]
    print(f"  implied major-element fraction: "
          + ("every class inside "
             f"{a['criteria']['implied_major_element_fraction']['range']}"
             if not outside else "OUTSIDE: " + "; ".join(outside)))

    b = result["by_drainage_fate"]
    print(f"\nRANK DIAGNOSTIC")
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
