#!/usr/bin/env python3
"""The land column property contract: one hydraulic and thermal description,
and the comparison of every consumer against it.

Worldbuilding. Vesper is an invented planet; everything below is about the
simulation of it -- a declared column of soil, the states and flow properties
it is described by, and what each of the three consumers currently derives for
itself instead.

    python pedology/scripts/land_column_properties.py            # status, exit 0
    python pedology/scripts/land_column_properties.py --strict   # exit 1 while any property is undeclared

This module is the enforcement for `pedology/config/land_column_properties.yaml`.
It does four things, and none of them runs a model:

1. **The declaration.** Units, geometry, ordering of the retention states, the
   potential convention and the uncertainty cases are checked for internal
   consistency, and the check is then run against declarations broken in named
   ways so it is falsifiable rather than merely quiet.

2. **The three-way capacity comparison**, on the checked-in soil map. Pedology's
   `awc`, ExoPlaSim's `dwmax` and LPJ-GUESS's own Cosby derivation with the
   Vesper regolith and bedrock scaling, on the same cells, in the same units.
   The audit's finding 1 is that these are not the same capacity; this is the
   number.

3. **The Cosby inversion's air-entry clamp.** `soilinput.cpp` clamps both
   matric heads at the air-entry value, because Cosby's equation holds only
   below air entry and unclamped it returns a field capacity above saturation.
   The texture region where the clamp binds is derived analytically from the
   regression coefficients and the soil map is measured against it, so the
   answer is a MARGIN rather than a count of zero.

4. **The gravity correction to field capacity.** Field capacity is a drainage
   equilibrium, so its matric pressure scales with gravity and this planet's is
   not Earth's. The wilting point is a plant pressure and does not move. The
   correction is therefore one-directional, has no free parameter, and neither
   consumer applies it.

The report goes to `pedology/analysis/land_column_properties_report.json`.
The contract is `pedology/notes/land-column-property-contract.md`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

from _paths import ANALYSIS, CONFIG, PROJECT_ROOT  # noqa: F401  (adds lib/ to sys.path)

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
DECLARATION = COMPONENT_ROOT / "config" / "land_column_properties.yaml"
REPORT = ANALYSIS / "land_column_properties_report.json"

UNDECLARED = "undeclared"

# ---------------------------------------------------------------------------
# LPJ-GUESS's own derivation, transcribed from the vendored source so the
# comparison runs against what the model computes rather than what a note says
# it computes. Every coefficient below is at `vendor/lpj-guess/modules/
# soilinput.cpp:360` in `SoilInput::get_mineral`, from Cosby et al. (1984)
# Table 4 and Equation 1.
#
# THIS IS A TRANSCRIPTION AND NOT A REIMPLEMENTATION. It exists to be compared
# against pedology's derivation on the same cells; it is not a second copy of
# the model's hydrology and nothing consumes it.
# ---------------------------------------------------------------------------

# `soilinput.cpp` works in RECIPROCAL suction: it sets Psi_s = 10^(-logPsi_s)
# where Cosby's regression gives log10 of the air-entry suction in cm of water,
# so a stored value of 10^-x is a suction of 10^x cm. Field capacity at 10^-2 is
# therefore 100 cm of water and the wilting point at 10^-4.2 is 15849 cm. The
# unit cancels out of the derivation, which enters only as (Psi_fc / Psi_s), so
# the states below are unit-free in the exponent. What they are NOT is free of
# gravity: 100 cm of water is a pressure only once a g is chosen, which is what
# `gravity_field_capacity_shift` below is about.
COSBY_PSI_FIELD_CAPACITY = 10.0 ** -2.0
COSBY_PSI_WILTING = 10.0 ** -4.2

# LPJ-GUESS's physical profile: `framework/guess.h` SOILDEPTH_UPPER = 500 mm
# over NSOILLAYER_UPPER = 5 layers, SOILDEPTH_LOWER = 1000 mm over the
# remaining 10 of NSOILLAYER = 15.
LPJ_SOILDEPTH_UPPER_MM = 500.0
LPJ_SOILDEPTH_LOWER_MM = 1000.0
LPJ_NSOILLAYER_UPPER = 5
LPJ_NSOILLAYER = 15


def load(path: Path = DECLARATION) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def cosby_parameters(sand, clay):
    """Cosby et al. (1984) Table 4 regressions, exactly as `get_mineral` uses them."""
    silt = 1.0 - sand - clay
    b = 3.10 + 15.7 * clay - 0.3 * sand
    log_psi_s = 1.54 - 0.95 * sand + 0.63 * silt
    theta_s = 0.01 * (50.5 - 14.2 * sand - 3.7 * clay)
    psi_s = 10.0 ** (-log_psi_s)
    return b, psi_s, theta_s


def cosby_states(sand, clay, psi_field_capacity=COSBY_PSI_FIELD_CAPACITY,
                 clamp=True):
    """Saturation, field capacity and wilting point, volumetric.

    CLAMPED AT AIR ENTRY, because the model is. Cosby's equation 1 holds only
    below the air-entry value; at and above it the pore space is full and
    theta is theta_s. `get_mineral` takes both heads through `min(..., Psi_s)`
    for that reason, so this transcription does too -- WORLD-NGA10 put the
    clamp in the code, and a transcription that kept the unclamped form would
    be comparing against a model that no longer exists.

    `clamp=False` is what `air_entry_clamp_region` uses to measure how far the
    soil map sits from the texture corner where the clamp starts binding. That
    margin is a property of the map and is worth carrying whether or not the
    code is right, because it is what says how close the nearest cell is.
    """
    b, psi_s, theta_s = cosby_parameters(sand, clay)
    psi_fc = np.minimum(psi_field_capacity, psi_s) if clamp else psi_field_capacity
    psi_wp = np.minimum(COSBY_PSI_WILTING, psi_s) if clamp else COSBY_PSI_WILTING
    theta_fc = theta_s * (psi_fc / psi_s) ** (1.0 / b)
    theta_wp = theta_s * (psi_wp / psi_s) ** (1.0 / b)
    return theta_s, theta_fc, theta_wp


def lpj_capacity_mm(sand, clay, regolith_depth_m, bedrock_fraction):
    """LPJ-GUESS's plant-available capacity over its whole profile, mm.

    Cosby field capacity minus wilting point, applied to the fifteen physical
    layers, then scaled layer by layer by the share above the regolith contact
    plus the bedrock fraction below it -- which is `vesperinput.cpp`'s
    `apply_regolith_scaling`, transcribed. The bedrock fraction is capped at 1
    there and is capped here.
    """
    _, theta_fc, theta_wp = cosby_states(sand, clay)
    available = theta_fc - theta_wp
    bedrock = np.minimum(np.asarray(bedrock_fraction, dtype=float), 1.0)
    depth_mm = np.asarray(regolith_depth_m, dtype=float) * 1000.0

    upper_layer_mm = LPJ_SOILDEPTH_UPPER_MM / LPJ_NSOILLAYER_UPPER
    lower_layer_mm = LPJ_SOILDEPTH_LOWER_MM / (LPJ_NSOILLAYER - LPJ_NSOILLAYER_UPPER)

    total = np.zeros_like(available, dtype=float)
    top_mm = np.zeros_like(available, dtype=float)
    for layer in range(LPJ_NSOILLAYER):
        thickness = upper_layer_mm if layer < LPJ_NSOILLAYER_UPPER else lower_layer_mm
        share = np.clip((depth_mm - top_mm) / thickness, 0.0, 1.0)
        usable = share + (1.0 - share) * bedrock
        total += available * thickness * usable
        top_mm = top_mm + thickness
    return total


# ---------------------------------------------------------------------------
# 1. The declaration
# ---------------------------------------------------------------------------

def check_declaration(decl: dict) -> list[str]:
    """Internal inconsistencies in the contract. Empty means it hangs together."""
    bad: list[str] = []

    geom = decl["geometry"]
    reach = geom["physical_layer_thickness_m"] * geom["physical_layer_count"]
    if abs(reach - geom["column_base_m"]) > 1e-9:
        bad.append(f"geometry: {geom['physical_layer_count']} layers of "
                   f"{geom['physical_layer_thickness_m']} m reach {reach} m, and "
                   f"column_base_m is {geom['column_base_m']}")
    if geom["rootable_base_m"] > geom["column_base_m"]:
        bad.append("geometry: the rootable base is below the column base. A root "
                   "cannot draw from water the column does not carry")
    if abs(geom["vadose_base_m"] - decl["aquifer_boundary"]["contact_depth_m"]) > 1e-9:
        bad.append("geometry: vadose_base_m and aquifer_boundary.contact_depth_m "
                   "are two names for one contact and they disagree")

    states = decl["states"]
    order = [s.strip() for s in states["ordering"].split("<=")]
    for name in order:
        if name not in states:
            bad.append(f"states: ordering names {name!r}, which is not a declared state")
    for name, body in states.items():
        if name == "ordering":
            continue
        if "gravity_dependent" not in body and "definition" in body:
            bad.append(f"state {name}: does not say whether it depends on gravity. "
                       "That is the one question this contract exists to answer "
                       "the same way twice")

    pot = decl["potential_convention"]
    if pot["wilting_point_pa"] >= 0:
        bad.append("potential_convention: the wilting point is a suction and is "
                   "negative in this contract's sign convention")
    if pot["field_capacity_drainage_length_m"] <= 0:
        bad.append("potential_convention: the drainage length is a length and is positive")
    # Cosby's own field-capacity suction IS rho_w * g_earth * 1 m, and that
    # coincidence is what licenses stating field capacity as a drainage length
    # at all. If the declared length stops reproducing it, the derivation below
    # has silently changed meaning.
    earth_fc_pa = (pot["water_density_kg_m3"] * pot["earth_gravity_m_s2"]
                   * pot["field_capacity_drainage_length_m"])
    if abs(earth_fc_pa - 9806.65) > 1.0:
        bad.append(f"potential_convention: the drainage length gives an Earth field "
                   f"capacity suction of {earth_fc_pa:.1f} Pa, and Cosby's 10^2 cm "
                   "of water is 9806.65 Pa. The two must agree or the drainage "
                   "length is not what Cosby measured at")

    closure = decl["retention_closure"]
    if closure.get("family") not in ("brooks_corey", "van_genuchten"):
        bad.append(f"retention_closure: family {closure.get('family')!r} is not one "
                   "of the two implemented families. A closure with no family stated "
                   "is the thing this block exists to prevent")
    if closure.get("name") in (None, UNDECLARED):
        bad.append("retention_closure: no name. A consumer has to be able to say "
                   "which curve it is reading")

    unc = decl["uncertainty"]
    if unc.get("applied") != "coherently across all layers and all properties in one case":
        bad.append("uncertainty: cases must be applied coherently. Drawing each "
                   "layer independently averages the spread away over 15 layers "
                   "and reports a false precision")
    if len(unc.get("cases", [])) < 2:
        bad.append("uncertainty: fewer than two cases is not a bracket")

    aq = decl["aquifer_boundary"]
    if not aq.get("scaled_separately"):
        bad.append("aquifer_boundary: transmissivity must be scaled separately. A "
                   "vadose unsaturated conductivity is not an aquifer permeability")

    for name, body in decl["units"].items():
        if not body:
            bad.append(f"units: {name} has no unit declared")
    return bad


def undeclared_properties(decl: dict) -> list[dict]:
    """Every property still carrying the sentinel, with the issue that owns it.

    ANY leaf equal to `undeclared` counts, not a fixed list of key names. The
    first version of this looked for `value` and `fraction_from` only, and went
    quiet the moment a property was refined enough to name its own field: the
    frozen-pore impedance gained a declared FORM with an undeclared EXPONENT and
    silently left the register while still being undeclared. A sentinel that
    only fires on the key names someone thought of is not a sentinel.
    """
    found: list[dict] = []

    def owner_of(chain):
        for node in reversed(chain):
            if isinstance(node, dict) and node.get("owner"):
                return node["owner"]
        return "unowned"

    def walk(node, path, chain):
        if isinstance(node, dict):
            for key, value in node.items():
                if value == UNDECLARED:
                    found.append({"property": ".".join(path),
                                  "field": key,
                                  "owner": owner_of(chain + [node])})
                else:
                    walk(value, path + [key], chain + [node])
        elif isinstance(node, list):
            for item in node:
                walk(item, path, chain)

    for section in ("states", "flow", "thermal", "materials", "uncertainty"):
        walk(decl[section], [section], [])
    # One row per property, not one per undeclared field under it.
    seen, unique = set(), []
    for row in found:
        if row["property"] in seen:
            continue
        seen.add(row["property"])
        unique.append(row)
    return unique


# ---------------------------------------------------------------------------
# The declaration check, against declarations built wrong
# ---------------------------------------------------------------------------
#
# `check_declaration` returning nothing proves nothing on its own: a check that
# cannot fail and a contract that is correct look identical from there.

def _break_layer_reach(decl):
    decl["geometry"]["physical_layer_count"] = 12


def _break_root_below_column(decl):
    decl["geometry"]["rootable_base_m"] = 2.0


def _break_contact_disagrees(decl):
    decl["aquifer_boundary"]["contact_depth_m"] = 3.0


def _break_state_ordering_name(decl):
    decl["states"]["ordering"] = "residual <= turgor <= field_capacity <= saturation"


def _break_gravity_silence(decl):
    decl["states"]["field_capacity"].pop("gravity_dependent")


def _break_wilting_sign(decl):
    decl["potential_convention"]["wilting_point_pa"] = 1.5468e6


def _break_drainage_length(decl):
    decl["potential_convention"]["field_capacity_drainage_length_m"] = 0.33


def _break_closure_family(decl):
    decl["retention_closure"]["family"] = "empirical"


def _break_uncertainty_independence(decl):
    decl["uncertainty"]["applied"] = "per layer, independently"


def _break_shared_aquifer_property(decl):
    decl["aquifer_boundary"]["scaled_separately"] = False


DECLARATION_MUTATIONS = (
    ("layer_reach_disagrees", _break_layer_reach),
    ("root_below_column_base", _break_root_below_column),
    ("contact_depth_disagrees", _break_contact_disagrees),
    ("undeclared_state_in_ordering", _break_state_ordering_name),
    ("state_silent_on_gravity", _break_gravity_silence),
    ("wilting_point_sign", _break_wilting_sign),
    ("drainage_length_not_cosbys", _break_drainage_length),
    ("closure_family_unnamed", _break_closure_family),
    ("uncertainty_drawn_per_layer", _break_uncertainty_independence),
    ("aquifer_property_shared", _break_shared_aquifer_property),
)


def run_declaration_mutations(path: Path = DECLARATION) -> tuple[list[dict], list[str]]:
    rows, defects = [], []
    for name, mutate in DECLARATION_MUTATIONS:
        broken = load(path)
        mutate(broken)
        caught = check_declaration(broken)
        rows.append({"mutation": name, "failures": len(caught),
                     "caught": bool(caught),
                     "first": caught[0] if caught else None})
        if not caught:
            defects.append(f"declaration mutation {name} was not caught. The check "
                           "passes a contract built wrong in that way")
    return rows, defects


# ---------------------------------------------------------------------------
# 2. The three-way capacity comparison
# ---------------------------------------------------------------------------

def read_soil_map(path: Path) -> dict[str, np.ndarray]:
    data = np.genfromtxt(path, names=True)
    return {name: np.asarray(data[name], dtype=float) for name in data.dtype.names}


def percentiles(values: np.ndarray) -> dict:
    q = np.percentile(values, [10, 50, 90])
    return {"p10": float(q[0]), "median": float(q[1]), "p90": float(q[2]),
            "mean": float(values.mean())}


def compare_consumers(decl: dict, soil: dict[str, np.ndarray]) -> dict:
    """Pedology's capacity against LPJ-GUESS's, on the same cells.

    ExoPlaSim's is pedology's divided by 1000 and installed as a bucket depth,
    so it is the same field in metres and is reported as such rather than as a
    third estimate: `build_surface_soil_water.py` reads the `awc` column and
    converts. What differs is not the number but what it is asked to mean --
    a plant-available capacity standing in for an entire land column.
    """
    pedology_mm = soil["awc"]
    lpj_mm = lpj_capacity_mm(soil["sand"], soil["clay"],
                             soil["depth"], soil["bedrockfrac"])
    ratio = lpj_mm / np.maximum(pedology_mm, 1e-12)
    difference = np.abs(lpj_mm - pedology_mm)
    correlation = float(np.corrcoef(pedology_mm, lpj_mm)[0, 1])
    return {
        "cells": int(pedology_mm.size),
        "pedology_awc_mm": percentiles(pedology_mm),
        "lpj_guess_derived_mm": percentiles(lpj_mm),
        "exoplasim_dwmax_m": percentiles(pedology_mm / 1000.0),
        "exoplasim_source": "the pedology awc column, converted mm to m; not a "
                            "third derivation",
        "lpj_over_pedology_ratio": percentiles(ratio),
        "absolute_difference_mm": percentiles(difference),
        "pearson_correlation": correlation,
        "reading": "The two fields correlate because both start from related "
                   "texture and depth. They are not the same capacity: one is a "
                   "declared endmember mixture with additive organic and "
                   "allophane terms, the other is a Cosby retention inversion "
                   "with no andic term and no bulk density.",
        "soil_carbon_is_zero": bool(np.all(soil["soilc"] == 0.0)),
        "soil_carbon_note": "the checked-in map carries zero soil carbon "
                            "everywhere, so the organic mixing paths of both "
                            "derivations are inactive and the divergence "
                            "reported here is the MINERAL one alone",
    }


# ---------------------------------------------------------------------------
# 3. The Cosby inversion's air-entry clamp, and the map's margin from it
# ---------------------------------------------------------------------------

def air_entry_clamp_region(soil: dict[str, np.ndarray]) -> dict:
    """Where the Cosby inversion's air-entry clamp binds, and how far the map is.

    Unclamped, `theta_fc = theta_s * (psi_fc / psi_s) ** (1 / b)` exceeds
    `theta_s` when `psi_fc > psi_s`, that is when `log10(psi_s^-1) > 2`.
    Substituting the regression for that log and silt = 1 - sand - clay gives

        1.58 * sand + 0.63 * clay < 0.17

    which is a straight line in the texture simplex and needs no search.

    `get_mineral` now clamps both heads at `Psi_s`, so a cell inside that
    region no longer returns a field capacity above saturation; it returns
    saturation, which is the physically right answer and is a DIFFERENT number
    from the one the retention curve would give. The map is still measured
    against the line as a MARGIN, because a count of zero says nothing about
    how close the nearest cell is, and because a map that moved inside would
    have its field capacity set by the clamp rather than by its texture.
    """
    sand, clay = soil["sand"], soil["clay"]
    theta_s, theta_fc, theta_wp = cosby_states(sand, clay)
    _, theta_fc_raw, theta_wp_raw = cosby_states(sand, clay, clamp=False)
    criterion = 1.58 * sand + 0.63 * clay
    inside = criterion < 0.17
    headroom = theta_s - theta_fc
    return {
        "criterion": "1.58 * sand + 0.63 * clay < 0.17 is where the air-entry "
                     "clamp binds; unclamped the inversion returns a field "
                     "capacity above saturation there",
        "criterion_bound": 0.17,
        "cells_inside": int(inside.sum()),
        "closest_cell_criterion": float(criterion.min()),
        "minimum_saturation_headroom_volumetric": float(headroom.min()),
        "any_field_capacity_above_saturation": bool((theta_fc > theta_s).any()),
        "any_wilting_point_above_field_capacity": bool((theta_wp > theta_fc).any()),
        "clamp_changes_any_cell": bool(
            np.any(theta_fc_raw != theta_fc) or np.any(theta_wp_raw != theta_wp)),
        "verdict": "the CODE keeps this right, and the map has never needed to. "
                   "WORLD-NGA10 clamped both heads at Psi_s in `get_mineral`, "
                   "so the inversion cannot return a field capacity above "
                   "saturation for any texture. No cell of this map is inside "
                   "the region, so the clamp changes no current soil property; "
                   "the margin is what says how far that is from being true by "
                   "accident.",
        "owner": "world-nga10",
    }


# ---------------------------------------------------------------------------
# 4. The gravity correction to field capacity
# ---------------------------------------------------------------------------

def gravity_field_capacity_shift(decl: dict, soil: dict[str, np.ndarray],
                                 gravity_m_s2: float) -> dict:
    """What this planet's gravity does to field capacity, and only to it.

    Field capacity is a drainage equilibrium over a stated length, so its
    matric suction is `rho_w * g * L` and scales with gravity. The wilting
    point is a plant pressure and does not. Saturation is pore geometry and
    does not. So the available capacity falls, by an amount set entirely by the
    retention exponent, with no free parameter and no fitted term.

    The size of the effect is small and its DIRECTION is not in doubt, which is
    the case for reporting it rather than rounding it away: it moves every cell
    the same way, so it does not average out over the map the way a scatter
    would.
    """
    pot = decl["potential_convention"]
    rho = pot["water_density_kg_m3"]
    length = pot["field_capacity_drainage_length_m"]
    earth_pa = rho * pot["earth_gravity_m_s2"] * length
    vesper_pa = rho * gravity_m_s2 * length

    sand, clay = soil["sand"], soil["clay"]
    b, _, theta_s = cosby_parameters(sand, clay)
    # theta(psi) ~ psi ** (-1/b) in magnitude, so the ratio of the two field
    # capacities is the pressure ratio to the power -1/b and every other factor
    # cancels. Nothing about Cosby's Psi_s survives into this number.
    ratio = (vesper_pa / earth_pa) ** (-1.0 / b)

    theta_s_e, theta_fc_e, theta_wp_e = cosby_states(sand, clay)
    theta_fc_v = theta_fc_e * ratio
    available_e = theta_fc_e - theta_wp_e
    available_v = theta_fc_v - theta_wp_e
    available_ratio = available_v / np.maximum(available_e, 1e-12)

    return {
        "gravity_m_s2": gravity_m_s2,
        "earth_field_capacity_suction_pa": earth_pa,
        "vesper_field_capacity_suction_pa": vesper_pa,
        "suction_ratio": vesper_pa / earth_pa,
        "field_capacity_ratio": percentiles(ratio),
        "available_capacity_ratio": percentiles(available_ratio),
        "direction": "down, on every cell. Stronger gravity drains a column to a "
                     "drier state at equilibrium; the wilting point is a plant "
                     "pressure and does not move, so plant-available water falls.",
        "applied_by_any_consumer": False,
        "applied_note": "pedology's declared endmember capacities are Earth field "
                        "observations and LPJ-GUESS's Cosby inversion evaluates at "
                        "the Earth suction, so neither carries this. It is a "
                        "correction this contract makes available and no consumer "
                        "has taken.",
        "bracket_note": "an upper bound on the size of the effect, not a "
                        "correction factor to apply as it stands: it holds the "
                        "Clapp-Hornberger exponent fixed while gravity moves, and "
                        "the within-texture-class variance Cosby reports is larger "
                        "than the shift.",
    }


# ---------------------------------------------------------------------------
# 5. What unifying the two derivations would move
# ---------------------------------------------------------------------------

def unification_cost(decl: dict, soil: dict[str, np.ndarray],
                     gravity_m_s2: float) -> dict:
    """What each consumer's capacity becomes if one derivation is made central.

    WORLD-OF6N is the removal this contract was built to make attributable, and
    it needs a decision that is about the world rather than about the code: WHAT
    THE CENTRAL CASE IS. This costs the candidates so that decision is a
    one-liner rather than an investigation, and it runs no model.

    Three candidates, and they are not symmetric. Two of them are readings of
    the contract as it already stands, because the contract already NAMES a
    retention closure and already declares field capacity as a drainage
    equilibrium at `rho_w * g * L`; taken literally that determines the states
    from texture and the only open part is which gravity the suction is
    evaluated at. The third would replace the declared closure with pedology's
    endmember mixture, which is not a retention curve at all and cannot be
    expressed in the declared closure -- so it is a change to the DECLARATION
    and not a change to a consumer.
    """
    pedology_mm = soil["awc"]
    lpj_mm = lpj_capacity_mm(soil["sand"], soil["clay"],
                             soil["depth"], soil["bedrockfrac"])
    # The whole per-layer scaling is linear in the available capacity, so the
    # gravity shift's per-cell ratio carries through the profile exactly.
    shift = gravity_field_capacity_shift(decl, soil, gravity_m_s2)
    pot = decl["potential_convention"]
    b, _, _ = cosby_parameters(soil["sand"], soil["clay"])
    fc_ratio = ((gravity_m_s2 / pot["earth_gravity_m_s2"]) ** (-1.0 / b))
    _, theta_fc_e, theta_wp_e = cosby_states(soil["sand"], soil["clay"])
    available_ratio = ((theta_fc_e * fc_ratio - theta_wp_e)
                       / np.maximum(theta_fc_e - theta_wp_e, 1e-12))
    lpj_vesper_mm = lpj_mm * available_ratio

    def moves(new, old):
        ratio = new / np.maximum(old, 1e-12)
        return {"median_mm_before": float(np.median(old)),
                "median_mm_after": float(np.median(new)),
                "land_mean_mm_before": float(old.mean()),
                "land_mean_mm_after": float(new.mean()),
                "ratio": percentiles(ratio)}

    candidates = {
        "closure_at_earth_suction": {
            "what": "the contract's declared closure -- Clapp-Hornberger on "
                    "Cosby Table 4 -- evaluated at Cosby's own field-capacity "
                    "suction, which is the pressure LPJ-GUESS already uses",
            "is_a_reading_of_the_contract": True,
            "central_case_mm": percentiles(lpj_mm),
            "exoplasim_dwmax": moves(lpj_mm, pedology_mm),
            "lpj_guess_capacity": moves(lpj_mm, lpj_mm),
            "deletes": "pedology's endmember volumetric capacities, its "
                       "additive organic term and its allophane term, from "
                       "build_soil.py's awc column",
        },
        "closure_at_vesper_suction": {
            "what": "the same closure at the suction the contract's own "
                    "potential convention declares, rho_w * g * L at this "
                    "planet's gravity. This is the contract read literally: "
                    "the closure is named and field capacity is declared as a "
                    "drainage equilibrium over a declared length",
            "is_a_reading_of_the_contract": True,
            "central_case_mm": percentiles(lpj_vesper_mm),
            "exoplasim_dwmax": moves(lpj_vesper_mm, pedology_mm),
            "lpj_guess_capacity": moves(lpj_vesper_mm, lpj_mm),
            "deletes": "the same, and additionally moves LPJ-GUESS off the "
                       "Earth suction it evaluates at today",
        },
        "pedology_endmember_mixture": {
            "what": "pedology's declared endmember volumetric capacities with "
                    "its additive organic and allophane terms. NOT a retention "
                    "curve: it carries no saturation, no air entry and no "
                    "exponent, so adopting it means replacing the contract's "
                    "named closure rather than choosing between consumers",
            "is_a_reading_of_the_contract": False,
            "central_case_mm": percentiles(pedology_mm),
            "exoplasim_dwmax": moves(pedology_mm, pedology_mm),
            "lpj_guess_capacity": moves(pedology_mm, lpj_mm),
            "deletes": "soilinput.cpp's get_mineral pedotransfer and "
                       "vesperinput.cpp's regolith rescaling, and leaves the "
                       "contract with a capacity and no retention curve, so "
                       "saturation and matric potential would have to come "
                       "from somewhere else",
        },
    }

    andic_mm = (float(_pedogenesis()["andisol"]["volumetric_capacity_allophane"])
                * soil["andic"] * soil["depth"] * 1000.0)
    gap_mm = np.abs(lpj_mm - pedology_mm)
    return {
        "candidates": candidates,
        "andic_term_mm": percentiles(andic_mm),
        "andic_cells": int((soil["andic"] > 0.0).sum()),
        "andic_share_of_pedology_awc": percentiles(
            andic_mm / np.maximum(pedology_mm, 1e-12)),
        "andic_is_not_the_disagreement": (
            "the allophane term LPJ-GUESS has no equivalent of contributes a "
            "median of nothing and a p90 of a few per cent of pedology's "
            "capacity, against a median absolute gap between the two "
            "derivations of tens of millimetres. What separates them is that "
            "the declared endmember capacities and the Cosby retention "
            "inversion are different numbers for the same texture, not that "
            "one carries a material the other does not"),
        "median_absolute_gap_mm": float(np.median(gap_mm)),
        "not_measurable_here": (
            "the change in land runoff ratio and in E over R. Both need a "
            "paired baseline, pedology/scripts/probe_runoff_response.py needs "
            "a climatology to back-solve potential evaporation from, and "
            "config/planet.yaml declares baseline_climatology null. On this "
            "build the capacity change can be costed and its climate "
            "consequence cannot"),
        "not_execution_verified": (
            "nothing on the LPJ-GUESS side. It does not build on this tree "
            "until a baseline climatology exists, so every LPJ number here is "
            "a transcription of the vendored source evaluated offline"),
        "owner": "world-of6n",
    }


def _pedogenesis() -> dict:
    return yaml.safe_load((COMPONENT_ROOT / "config" / "pedogenesis.yaml")
                          .read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------

def build_report(decl: dict, soil_map: Path, gravity_m_s2: float) -> dict:
    soil = read_soil_map(soil_map)
    failures = check_declaration(decl)
    mutations, defects = run_declaration_mutations()
    return {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "declaration": "pedology/config/land_column_properties.yaml",
        "declaration_sha256": hashlib.sha256(DECLARATION.read_bytes()).hexdigest(),
        "version": decl["version"],
        "soil_map": str(soil_map.resolve().relative_to(PROJECT_ROOT)),
        "soil_map_sha256": hashlib.sha256(soil_map.read_bytes()).hexdigest(),
        "units": decl["units"],
        "geometry": decl["geometry"],
        "retention_closure": decl["retention_closure"],
        "declaration_failures": failures,
        "declaration_mutations": mutations,
        "checker_defects": defects,
        "consumer_comparison": compare_consumers(decl, soil),
        "cosby_air_entry_clamp": air_entry_clamp_region(soil),
        "gravity_field_capacity_shift": gravity_field_capacity_shift(
            decl, soil, gravity_m_s2),
        "unification_cost": unification_cost(decl, soil, gravity_m_s2),
        "undeclared_properties": undeclared_properties(decl),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--soil-map", type=Path, default=None,
                        help="defaults to the configured build's soil map")
    parser.add_argument("--strict", action="store_true",
                        help="exit 1 while any contract property is still undeclared")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    gravity = float(config["planet"]["gravity_m_s2"])
    if args.soil_map is None:
        sys.path.insert(0, str(PROJECT_ROOT / "lib"))
        import builds
        args.soil_map = builds.soilmap(config)
    if not args.soil_map.is_file():
        raise SystemExit(f"{args.soil_map} does not exist. Run "
                         "pedology/scripts/build_soil.py.")

    decl = load()
    report = build_report(decl, args.soil_map, gravity)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2, sort_keys=False) + "\n",
                      encoding="utf-8")

    hard = report["declaration_failures"] + report["checker_defects"]
    if not args.quiet:
        geom = decl["geometry"]
        print(f"land column properties v{report['version']}: "
              f"{geom['physical_layer_count']} physical layers to "
              f"{geom['column_base_m']} m, rootable to {geom['rootable_base_m']} m, "
              f"closure {decl['retention_closure']['name']} "
              f"({decl['retention_closure']['family']})")
        caught = sum(1 for row in report["declaration_mutations"] if row["caught"])
        print(f"  declaration mutations   {caught} of "
              f"{len(report['declaration_mutations'])} contracts built wrong were caught")

        cmp_ = report["consumer_comparison"]
        print(f"\n  capacity on {cmp_['cells']} cells of "
              f"{report['soil_map']}")
        print(f"    pedology awc          {cmp_['pedology_awc_mm']['median']:7.1f} mm "
              f"median, {cmp_['pedology_awc_mm']['p10']:.1f} to "
              f"{cmp_['pedology_awc_mm']['p90']:.1f} p10-p90")
        print(f"    lpj-guess derived     {cmp_['lpj_guess_derived_mm']['median']:7.1f} mm "
              f"median, {cmp_['lpj_guess_derived_mm']['p10']:.1f} to "
              f"{cmp_['lpj_guess_derived_mm']['p90']:.1f} p10-p90")
        print(f"    ratio lpj/pedology    {cmp_['lpj_over_pedology_ratio']['median']:7.2f} "
              f"median, {cmp_['lpj_over_pedology_ratio']['p10']:.2f} to "
              f"{cmp_['lpj_over_pedology_ratio']['p90']:.2f} p10-p90")
        print(f"    absolute difference   "
              f"{cmp_['absolute_difference_mm']['median']:7.1f} mm median")
        print(f"    pearson               {cmp_['pearson_correlation']:7.3f}")

        clamp = report["cosby_air_entry_clamp"]
        print(f"\n  cosby air-entry clamp   {clamp['cells_inside']} cells inside; "
              f"closest cell at {clamp['closest_cell_criterion']:.3f} against a "
              f"bound of {clamp['criterion_bound']}")
        print(f"    saturation headroom   "
              f"{clamp['minimum_saturation_headroom_volumetric']:.4f} volumetric at "
              "the closest cell")

        unify = report["unification_cost"]
        print("\n  what one central case would move, world-of6n")
        for name, cand in unify["candidates"].items():
            reading = "contract as written" if cand["is_a_reading_of_the_contract"] \
                else "replaces the declared closure"
            print(f"    {name:28s} ({reading})")
            for who in ("exoplasim_dwmax", "lpj_guess_capacity"):
                move = cand[who]
                print(f"      {who:20s} "
                      f"{move['median_mm_before']:7.1f} -> "
                      f"{move['median_mm_after']:7.1f} mm median, "
                      f"x{move['ratio']['median']:.3f} "
                      f"({move['ratio']['p10']:.3f} to {move['ratio']['p90']:.3f})")

        grav = report["gravity_field_capacity_shift"]
        print(f"\n  gravity shift at {grav['gravity_m_s2']} m/s2: field capacity "
              f"suction {grav['suction_ratio']:.3f}x Earth's")
        print(f"    field capacity        "
              f"{grav['field_capacity_ratio']['median']:.4f} median of Earth's, "
              f"{grav['field_capacity_ratio']['p10']:.4f} to "
              f"{grav['field_capacity_ratio']['p90']:.4f}")
        print(f"    available capacity    "
              f"{grav['available_capacity_ratio']['median']:.4f} median of Earth's, "
              f"{grav['available_capacity_ratio']['p10']:.4f} to "
              f"{grav['available_capacity_ratio']['p90']:.4f}")

        print(f"\n  {len(report['undeclared_properties'])} properties undeclared; "
              "the contract is defined and is not complete")
        for row in report["undeclared_properties"]:
            print(f"    {row['property'] + '.' + row['field']:<48} {row['owner']}")
        for problem in hard:
            print(f"  FAIL {problem}")
        print(f"\nwrote {REPORT.relative_to(PROJECT_ROOT)}")

    if hard:
        return 1
    if args.strict and report["undeclared_properties"]:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
