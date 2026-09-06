#!/usr/bin/env python3
"""The land column property contract: one hydraulic description of the simulated
land column, emitted for every consumer to read.

Worldbuilding. Vesper is an invented planet; everything below is about the
simulation of it -- a declared column of soil and the states and flow
properties it is described by.

    python pedology/scripts/land_column_properties.py            # check, emit, report
    python pedology/scripts/land_column_properties.py --strict   # exit 1 while any property is undeclared
    python pedology/scripts/land_column_properties.py --no-emit  # check and report only
    python pedology/scripts/land_column_properties.py --update-declaration
                                        # write the derived saturation mapping
                                        # back into the contract, then report
                                        # every restatement of it that is stale

This module is the enforcement for `pedology/config/land_column_properties.yaml`
AND the one implementation of the states it declares. It does five things, and
none of them runs a model:

1. **The declaration.** Units, geometry, ordering of the retention states, the
   potential convention and the uncertainty cases are checked for internal
   consistency, and the check is then run against declarations broken in named
   ways so it is falsifiable rather than merely quiet.

2. **The emission.** The adopted per-cell retention states and the per-layer
   weathered-bedrock fractions, written to the build's states file. ExoPlaSim
   installs the capacity column as `dwmax`; LPJ-GUESS takes the states and the
   fractions through its driver file. Neither derives them, and the derivation
   is here and nowhere else.

3. **The frame check.** Cosby's parameters are recorded as heads on Earth, and
   a head is a pressure only through the local gravity. Air entry is a
   capillary pressure and the wilting point is a plant pressure, so both are
   invariant and both of their heads scale together; field capacity is a
   drainage equilibrium and is the one state whose defining pressure carries
   gravity. The check computes the states in two consistent frames, which must
   agree exactly, and in the frame-mixed variant, which must not.

4. **The Cosby inversion's air-entry clamp.** Cosby's equation holds only below
   air entry and unclamped it returns a field capacity above saturation. The
   texture region where the clamp binds is derived analytically from the
   regression coefficients and the soil map is measured against it, so the
   answer is a MARGIN rather than a count of zero.

5. **What the adopted case moved**, against pedology's declared endmember
   mixture and against the suction the shipped LPJ-GUESS code evaluated at.
   That is the attribution the removal of the consumer-side pedotransfers
   needed.

The report goes to `pedology/analysis/land_column_properties_report.json` and
carries the emitted file's hash, so the emission's provenance is in it.
The contract is `pedology/notes/land-column-property-contract.md`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

from _paths import ANALYSIS, CONFIG, PROJECT_ROOT  # noqa: F401  (adds lib/ to sys.path)

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
DECLARATION = COMPONENT_ROOT / "config" / "land_column_properties.yaml"
REPORT = ANALYSIS / "land_column_properties_report.json"

# THE FILE THAT RESTATES THE SATURATION MAPPING IN COMPILED FORM. Read here,
# never written: `scan_restatements` says which copies of the declared pair have
# gone stale, and each copy has its own enforcer -- `run_exoplasim.py` refuses at
# import when landmod's compiled defaults disagree with this contract, and
# `scripts/check_consistency.py` refuses when the project config does.
LANDMOD_SOURCE = (PROJECT_ROOT / "vendor" / "exoplasim" / "exoplasim"
                  / "plasim" / "src" / "landmod.f90")

UNDECLARED = "undeclared"

# ---------------------------------------------------------------------------
# THE SUPERSEDED LPJ-GUESS DERIVATION, kept as a baseline and not as a model.
#
# `SoilInput::get_mineral` inverted Cosby et al. (1984) Table 4 and Equation 1
# for itself, at Cosby's own field-capacity and wilting suctions, and
# `VesperInput::apply_regolith_depth` then rescaled the profile. WORLD-OF6N
# removed both: the model reads `contract_states` and `layer_usable_fraction`
# above. What survives here is what those two computed, because the cost of the
# replacement is only legible against the number it replaced.
#
# Cosby's Table 4 regressions themselves are NOT superseded -- they are the
# contract's declared `parameter_source`, and `cosby_parameters` is the one
# place they are written. What is superseded is the SUCTION they were evaluated
# at, which was Earth's.
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


def cosby_parameter_spread(sand, clay):
    """Cosby Table 4's S.D. rows: the residual spread WITHIN a texture class.

    THE UNCERTAINTY OF THE ADOPTED CENTRAL CASE, and it comes from the same
    table and the same fit as the central case itself. Cosby et al. regressed
    the mean of each hydraulic parameter on texture, which is what
    `cosby_parameters` evaluates, and they regressed the STANDARD DEVIATION of
    each parameter within each of the eleven textural classes on the same
    texture. The second regression is what says how much of a soil a texture
    does not determine.

    Their coefficients are per PERCENT of sand, silt and clay and are converted
    to fractions here, once, the way every other Earth-recorded constant in this
    file is converted at the point it is declared:

        S.D. b        = 0.92 + 0.0492 %clay + 0.0144 %silt   R2 0.524, p 0.012
        S.D. log Psi_s= 0.72 - 0.0026 %silt + 0.0012 %clay   R2 0.096, p 0.355
        S.D. Theta_s  = 8.23 - 0.0805 %clay - 0.0070 %sand   R2 0.567, p 0.007

    Two limits on reading these, both from Cosby's own table. The second
    variable in each row is not significant, exactly as it is not in the mean
    rows the central case uses, and the `log Psi_s` row explains almost none of
    the between-class variation in its own spread -- so that spread is close to
    a constant 0.7 in log10 and is used as one number rather than as a texture
    response. What none of them carries is the CORRELATION between the three
    parameters' errors, which Cosby does not report, and that is why the cases
    built from these are an envelope and not a distribution.

    The same table carries `S.D. log Ks`, which nothing here reads: this
    contract declares no saturated conductivity, so the spread of one would be
    the spread of a number that does not exist. It arrives with `lshy-3`.
    """
    sand = np.asarray(sand, dtype=float)
    clay = np.asarray(clay, dtype=float)
    silt = 1.0 - sand - clay
    return {
        "b": 0.92 + 4.92 * clay + 1.44 * silt,
        "log_psi_s": 0.72 - 0.26 * silt + 0.12 * clay,
        "theta_s": 0.01 * (8.23 - 8.05 * clay - 0.70 * sand),
    }


def cosby_parameters(sand, clay, deviates=None):
    """Cosby et al. (1984) Table 4 regressions. The contract's parameter source.

    `psi_s` is returned in `get_mineral`'s old reciprocal convention -- a stored
    10^-x is a suction of 10^x cm of water -- because the two derivations below
    are both written against it. `air_entry_pressure_pa` is what turns it into
    the pressure the adopted closure works in.

    `deviates` is `(z_b, z_log_psi_s, z_theta_s)` in units of the within-class
    standard deviation `cosby_parameter_spread` gives, and `None` is the central
    case. It moves the PARAMETERS and not the states, so every state, every
    layer and every consumer of one case sees the same soil: that is what the
    contract means by applying a case coherently.
    """
    silt = 1.0 - sand - clay
    b = 3.10 + 15.7 * clay - 0.3 * sand
    log_psi_s = 1.54 - 0.95 * sand + 0.63 * silt
    theta_s = 0.01 * (50.5 - 14.2 * sand - 3.7 * clay)
    if deviates is not None:
        spread = cosby_parameter_spread(sand, clay)
        z_b, z_psi, z_theta = deviates
        b = b + z_b * spread["b"]
        log_psi_s = log_psi_s + z_psi * spread["log_psi_s"]
        theta_s = theta_s + z_theta * spread["theta_s"]
        # A case has to stay a soil. The exponent of a retention curve is
        # positive and a porosity is a fraction of a bulk volume; a case that
        # leaves either is not a wider bracket, it is a different closure.
        if np.any(b <= 0.0) or np.any(theta_s <= 0.0) or np.any(theta_s >= 1.0):
            raise SystemExit(
                f"uncertainty case {deviates} puts a cell outside the closure: "
                f"b down to {float(np.min(b)):.3f} and porosity in "
                f"[{float(np.min(theta_s)):.3f}, {float(np.max(theta_s)):.3f}]. "
                "The within-class spread is a spread of Cosby's parameters and "
                "cannot be applied where it leaves the parameter space.")
    psi_s = 10.0 ** (-log_psi_s)
    return b, psi_s, theta_s


def cosby_states(sand, clay, psi_field_capacity=COSBY_PSI_FIELD_CAPACITY,
                 clamp=True):
    """Saturation, field capacity and wilting point at COSBY'S OWN suction.

    THE SUPERSEDED DERIVATION, kept as the baseline the adopted case is costed
    against. `contract_states` is what the world's soil is described by; this
    evaluates the same closure at 10^2 cm of water, which is a pressure only at
    Earth's gravity, and is what `soilinput.cpp` computed for itself before it
    read the contract.

    CLAMPED AT AIR ENTRY. Cosby's equation 1 holds only below the air-entry
    value; at and above it the pore space is full and theta is theta_s. This is
    the shipped LPJ-GUESS derivation as it stood when WORLD-OF6N replaced it,
    kept so the report can say what the replacement moved, and the clamp is
    part of what it was.

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
    """The SHIPPED LPJ-GUESS capacity over its whole profile, mm.

    Cosby field capacity minus wilting point at Cosby's own suction, applied to
    the fifteen physical layers, then scaled layer by layer by the share above
    the regolith contact plus the bedrock fraction below it. Both halves were
    consumer-side derivations and both are gone; this is what they produced,
    kept so the adopted case has a baseline.
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
# THE ADOPTED DERIVATION. One implementation, in pressure, and the only place
# the simulated column's retention states are computed.
#
# The contract names its closure -- Clapp-Hornberger on Cosby Table 4 -- and
# declares field capacity as a drainage equilibrium at `rho_w * g * L`. Read
# literally that determines the states from texture, and this is that reading.
# WORLD-OF6N adopted it; WORLD-SLPA is the same change reached from the gravity
# side and is subsumed by it.
#
# WHY THIS WORKS IN PASCALS AND NOT IN HEADS. Cosby's parameters are recorded
# as heads in centimetres of water, and a head is a pressure only through the
# local gravity. Two of the three defining pressures are INVARIANT under a
# change of gravity and one is not:
#
#   the air-entry value is a CAPILLARY pressure, set by pore geometry and
#   surface tension. Invariant.
#   the wilting point is a PLANT pressure, the suction a root can generate.
#   Invariant.
#   field capacity is a drainage equilibrium against the weight of a column,
#   `rho_w * g * L`. It is the only one that carries gravity.
#
# Converting one of the two invariant heads and not the other mixes two frames
# and produces a wilting point that moves for a bookkeeping reason. Working in
# pressure makes that impossible to write: each Earth-recorded head is turned
# into a pressure ONCE, at Earth's gravity, at the point it is declared, and
# this world's gravity enters at exactly one place, the field-capacity
# equilibrium. `frame_consistency` below is the check, and it has something
# that can fail.
# ---------------------------------------------------------------------------

CM_WATER_PER_M = 100.0

# The guard `vesperinput.cpp` applies, kept identical to it because the
# weathered-bedrock rule is one rule and this is now the copy that computes it.
# LPJ-GUESS carries soil water as a fraction of a layer's capacity, so a layer
# scaled to exactly zero divides by zero and the NaN zeroes the gridcell's
# vegetation silently rather than failing.
NUMERICAL_FLOOR = 1.0e-4


def air_entry_pressure_pa(sand, clay, decl: dict, deviates=None):
    """Cosby's air-entry head as a pressure. Earth-recorded, and invariant.

    `log_psi_s` is log10 of the air-entry suction in centimetres of water as
    Cosby's Table 4 regression gives it, which is a head measured on Earth. It
    becomes a pressure through Earth's gravity, once, here.
    """
    pot = decl["potential_convention"]
    _, psi_s, _ = cosby_parameters(sand, clay, deviates)
    head_cm = 1.0 / psi_s            # the suction itself, cm of water on Earth
    return (pot["water_density_kg_m3"] * pot["earth_gravity_m_s2"]
            * head_cm / CM_WATER_PER_M)


def field_capacity_pressure_pa(decl: dict, gravity_m_s2: float) -> float:
    """The drainage equilibrium this world's column settles at.

    The ONE place gravity enters the retention states. `L` is the declared
    drainage length; the pressure that balances the weight of that column is
    `rho_w * g * L`, and `g` is this planet's.
    """
    pot = decl["potential_convention"]
    return (pot["water_density_kg_m3"] * gravity_m_s2
            * pot["field_capacity_drainage_length_m"])


def contract_states(sand, clay, decl: dict, gravity_m_s2: float, deviates=None):
    """The adopted volumetric states: saturation, field capacity, wilting point.

    Clapp-Hornberger, `theta(p) = theta_s * (p_ae / p) ** (1 / b)`, clamped at
    air entry because Cosby's equation holds only below it: at and above it the
    pore space is full and theta is theta_s. Returns `b` as well, because the
    exponent is a declared property of the closure and consumers that need a
    texture-dependent drainage rule take it from here rather than refitting one.

    `deviates` selects an uncertainty case; `None` is the central case and is
    what every consumer installs. See `cosby_parameters`.
    """
    pot = decl["potential_convention"]
    b, _, theta_s = cosby_parameters(sand, clay, deviates)
    p_ae = air_entry_pressure_pa(sand, clay, decl, deviates)
    p_fc = np.maximum(field_capacity_pressure_pa(decl, gravity_m_s2), p_ae)
    p_wilt = np.maximum(abs(float(pot["wilting_point_pa"])), p_ae)
    theta_fc = theta_s * (p_ae / p_fc) ** (1.0 / b)
    theta_wp = theta_s * (p_ae / p_wilt) ** (1.0 / b)
    return b, theta_s, theta_fc, theta_wp


def layer_usable_fraction(depth_m, bedrock_fraction, decl: dict):
    """Per physical layer, the share of a layer's capacity the column carries.

    The contract's `materials.weathered_bedrock` vertical rule, and the ONE
    implementation of it. A layer entirely above the regolith contact keeps all
    of its capacity; the part below the contact keeps `bedrock_fraction` of it,
    capped at the declared `capacity_ceiling` so sub-bedrock material can at
    most match the soil above.

    Returned as a fraction per layer rather than as a capacity, because the same
    fraction scales this layer's available water, its wilting point and its
    saturation: the material below the contact holds less of everything, and
    emitting one number instead of three is what stops them drifting apart.

    `NUMERICAL_FLOOR` is a numerical guard and not physics; it is declared with
    the constant.
    """
    geom = decl["geometry"]
    nlayer = int(geom["physical_layer_count"])
    thickness_mm = float(geom["physical_layer_thickness_m"]) * 1000.0
    ceiling = float(decl["materials"]["weathered_bedrock"]["capacity_ceiling"])

    depth_mm = np.atleast_1d(np.asarray(depth_m, dtype=float)) * 1000.0
    bedrock = np.atleast_1d(np.asarray(bedrock_fraction, dtype=float))
    bedrock = np.clip(bedrock, NUMERICAL_FLOOR, ceiling)

    usable = np.empty((depth_mm.size, nlayer), dtype=float)
    top_mm = 0.0
    for layer in range(nlayer):
        share = np.clip((depth_mm - top_mm) / thickness_mm, 0.0, 1.0)
        usable[:, layer] = share + (1.0 - share) * bedrock
        top_mm += thickness_mm
    return usable


def column_capacity_mm(theta_fc, theta_wp, usable, decl: dict):
    """Plant-available capacity over the whole physical column, mm.

    The states integrated over the declared layers with the weathered-bedrock
    rule applied. This is the capacity ExoPlaSim installs as `dwmax` and the
    same one LPJ-GUESS's layers sum to, which is what makes the two consumers
    agree on how much soil there is by construction rather than by coincidence.
    """
    geom = decl["geometry"]
    thickness_mm = float(geom["physical_layer_thickness_m"]) * 1000.0
    available = np.atleast_1d(np.asarray(theta_fc, dtype=float)
                              - np.asarray(theta_wp, dtype=float))
    return available * thickness_mm * usable.sum(axis=1)


def surface_layer_capacity_mm(theta_fc, theta_wp, usable, decl: dict):
    """(capacity, sub-wilting increment) of the surface layer, both in mm.

    The contract's `surface_layer` cut from AIR DRY rather than from the
    wilting point, which is the one place this column's capacity is not
    plant-available: a root cannot reach the water between air dry and the
    wilting point and a bare drying surface can. Returns the increment
    separately because that is exactly the difference between what LPJ-GUESS
    reads and what ExoPlaSim installs, and a difference nobody can name is the
    failure mode this contract exists to prevent.

    The layer's usable fraction is the FIRST physical layer's, because the
    surface layer is inside it: 0.02 m sits within the top 0.1 m whatever the
    regolith contact does.
    """
    surf = decl["surface_layer"]
    thickness_mm = float(surf["thickness_m"]) * 1000.0
    air_dry = float(surf["air_dry_water_content"])
    theta_fc = np.atleast_1d(np.asarray(theta_fc, dtype=float))
    theta_wp = np.atleast_1d(np.asarray(theta_wp, dtype=float))
    top = usable[:, 0]
    capacity = (theta_fc - air_dry) * thickness_mm * top
    increment = (theta_wp - air_dry) * thickness_mm * top
    return capacity, increment


def evaporable_capacity_mm(theta_fc, theta_wp, usable, decl: dict):
    """The capacity ExoPlaSim installs as `dwmax`, mm.

    `awc_mm` plus the surface layer's sub-wilting water, and nothing else. The
    identity `evaporable_mm - awc_mm == surface_layer increment` holds cell by
    cell by construction and is checked in `surface_layer_report`, which is the
    conservation statement across the split: the column gained exactly the
    water the surface layer can now reach and not a millimetre more.
    """
    awc = column_capacity_mm(theta_fc, theta_wp, usable, decl)
    _capacity, increment = surface_layer_capacity_mm(theta_fc, theta_wp, usable, decl)
    return awc + increment


def cascade_step(capacities, store, flux_mm):
    """`landcolumn.f90:column_step` transcribed, impermeable base, mm.

    A transcription and not a second scheme: the fill-and-pass order, the
    downward draw on a negative top layer, the base backing up and the floor at
    zero are the Fortran's, operation for operation. It is here so that the
    invariance the surface layer rests on can be DRIVEN rather than asserted.
    """
    w = list(store)
    w[0] += flux_mm
    for j in range(len(w) - 1):
        if w[j] < 0.0:
            w[j + 1] += w[j]
            w[j] = 0.0
    for j in range(len(w) - 1):
        excess = max(0.0, w[j] - capacities[j])
        w[j] = min(capacities[j], w[j])
        w[j + 1] += excess
    excess = max(0.0, w[-1] - capacities[-1])
    w[-1] = min(capacities[-1], w[-1])
    for j in range(len(w) - 2, -1, -1):
        w[j] += excess
        excess = max(0.0, w[j] - capacities[j])
        w[j] = min(capacities[j], w[j])
    return [max(v, 0.0) for v in w], excess


def check_split_invariance(seed: int = 20260826, steps: int = 4096) -> dict:
    """THE CHECK WITH A RIGHT ANSWER: splitting a layer moves no water.

    The surface layer is a repartition of the top of the climate column, and
    the claim it rests on is that the cascade's TOTAL is invariant under any
    repartition of a fixed total capacity: what enters is the surface flux,
    what leaves is base overflow, and neither depends on where the internal
    boundaries sit. Drive a two-layer column and a three-layer column whose
    first two capacities sum to the two-layer column's first, with the SAME
    flux sequence, and the total and the runoff must agree exactly.

    This is a conservation law and not a comparison: the answer is equality,
    and any difference at all is a defect. Reported as the largest absolute
    disagreement over the sequence, in mm, against a bound of zero plus the
    accumulated floating-point rounding of the sums, which is what
    `total_tolerance_mm` is.
    """
    rng = np.random.default_rng(seed)
    cap2 = [180.0, 300.0]
    cap3 = [3.0, 177.0, 300.0]
    assert abs(sum(cap2) - sum(cap3)) < 1e-12
    w2 = [90.0, 150.0]
    w3 = [3.0, 87.0, 150.0]
    worst_total = 0.0
    worst_runoff = 0.0
    for flux in rng.normal(0.0, 12.0, steps):
        w2, off2 = cascade_step(cap2, w2, float(flux))
        w3, off3 = cascade_step(cap3, w3, float(flux))
        worst_total = max(worst_total, abs(sum(w2) - sum(w3)))
        worst_runoff = max(worst_runoff, abs(off2 - off3))
    tolerance = 1.0e-9
    return {
        "what": "the column total and the runoff under a two-layer and a "
                "three-layer partition of the same total capacity, driven by "
                "one flux sequence through landcolumn.f90's own cascade",
        "steps": steps,
        "worst_total_disagreement_mm": float(worst_total),
        "worst_runoff_disagreement_mm": float(worst_runoff),
        "total_tolerance_mm": tolerance,
        "passes": bool(worst_total <= tolerance and worst_runoff <= tolerance),
        "why_it_matters": "the surface layer buys a store that empties in a day "
                          "instead of a season and costs the water budget "
                          "nothing. If this fails, it costs something, and the "
                          "cost is not a rounding error.",
    }


def frame_consistency(soil: dict[str, np.ndarray], decl: dict,
                      gravity_m_s2: float) -> dict:
    """Does the adopted derivation depend on which units the heads are carried in?

    IT MUST NOT, and this is the check with a right answer rather than a
    difference to eyeball. The same three states are computed three ways:

    - EARTH-EQUIVALENT FRAME. Every head stays in centimetres of water on
      Earth, and the field-capacity head is raised by `g / g_earth` so that it
      still names the pressure a column of the declared length exerts here.
    - VESPER FRAME. Every Earth-recorded head is divided by `g / g_earth` so it
      names the same pressure as a head of water on this world, and the
      field-capacity head is the declared drainage length in that unit.
    - FRAME-MIXED. The wilting head converted and the air-entry head left
      alone. This is not a frame at all, and it is here because a check that
      cannot fail proves nothing.

    The first two must agree exactly: `theta = theta_s * (h_ae / h) ** (1 / b)`
    depends on the two heads only through their ratio, and a frame change
    multiplies both by the same factor. The third must NOT agree, and the size
    of its disagreement is what a consumer would have shipped had it converted
    the plant pressure and forgotten the capillary one.
    """
    pot = decl["potential_convention"]
    sand, clay = soil["sand"], soil["clay"]
    b, _, theta_s = cosby_parameters(sand, clay)
    factor = gravity_m_s2 / pot["earth_gravity_m_s2"]

    head_ae_earth = 1.0 / cosby_parameters(sand, clay)[1]
    head_wilt_earth = (abs(float(pot["wilting_point_pa"]))
                       / (pot["water_density_kg_m3"] * pot["earth_gravity_m_s2"])
                       * CM_WATER_PER_M)
    # The declared drainage length as a head of water on Earth, which is what
    # Cosby's field-capacity suction is and why the length is declared at all.
    head_fc_earth = (float(pot["field_capacity_drainage_length_m"])
                     * CM_WATER_PER_M)

    def available_at(head_fc, head_wilt, head_ae):
        theta_fc = theta_s * (head_ae / np.maximum(head_fc, head_ae)) ** (1.0 / b)
        theta_wp = theta_s * (head_ae / np.maximum(head_wilt, head_ae)) ** (1.0 / b)
        return theta_fc - theta_wp

    shipped = available_at(head_fc_earth, head_wilt_earth, head_ae_earth)
    earth_frame = available_at(head_fc_earth * factor, head_wilt_earth,
                               head_ae_earth)
    vesper_frame = available_at(head_fc_earth, head_wilt_earth / factor,
                                head_ae_earth / factor)
    mixed = available_at(head_fc_earth, head_wilt_earth / factor, head_ae_earth)

    _, _, theta_fc_a, theta_wp_a = contract_states(sand, clay, decl, gravity_m_s2)
    adopted = theta_fc_a - theta_wp_a

    disagreement = float(np.max(np.abs(earth_frame - vesper_frame)))
    from_adopted = float(np.max(np.abs(adopted - earth_frame)))
    tolerance = 1e-12
    defects = []
    if disagreement > tolerance:
        defects.append(
            f"frame consistency: the Earth-equivalent and Vesper frames differ "
            f"by {disagreement:.3e} volumetric and must agree exactly. The "
            "states depend on the heads only through their ratio, so a frame "
            "change cancels; a difference means one head was converted and "
            "another was not")
    if from_adopted > tolerance:
        defects.append(
            f"frame consistency: the adopted pressure derivation differs from "
            f"the head derivation by {from_adopted:.3e} volumetric. They are "
            "the same closure written twice and must agree")
    if float(np.max(np.abs(mixed - earth_frame))) <= tolerance:
        defects.append(
            "frame consistency: the frame-mixed variant agrees with the "
            "consistent frames, so this check cannot tell them apart and "
            "proves nothing")

    def ratio(new):
        return percentiles(new / np.maximum(shipped, 1e-12))

    return {
        "gravity_ratio_to_earth": float(factor),
        "field_capacity_head_cm_earth_frame": float(head_fc_earth * factor),
        "available_ratio_to_shipped": {
            "earth_equivalent_frame": ratio(earth_frame),
            "vesper_frame": ratio(vesper_frame),
            "frame_mixed_wilting_only": ratio(mixed),
        },
        "frames_agree_max_abs_volumetric": disagreement,
        "adopted_matches_head_derivation_max_abs_volumetric": from_adopted,
        "frame_mixed_is_distinguishable": bool(
            float(np.max(np.abs(mixed - earth_frame))) > tolerance),
        "reading": "the two consistent frames are one number and the "
                   "frame-mixed variant is a different one. The wilting point "
                   "does not move under a change of gravity, because it and "
                   "the air-entry value are both invariant pressures recorded "
                   "as Earth heads and their ratio is what the closure reads. "
                   "Field capacity is the only state whose defining pressure "
                   "carries gravity.",
        "defects": defects,
        "owner": "world-slpa, subsumed into world-of6n",
    }


# ---------------------------------------------------------------------------
# THE DERIVED SCALARS, AND WHO WRITES THEM
#
# `thermal.saturation_mapping`'s four numbers are the median and the spread of a
# distribution this contract emits per cell, so they move whenever the soil map
# is rebuilt, which is every turn of loop A. Typed by hand they drift silently
# in one direction: the declaration goes on describing the previous build's soil
# while `check_thermal_against_states` is the only thing that knows. So they are
# WRITTEN from the same measurement that check makes, by `--update-declaration`,
# and the check is what proves the declaration is that measurement rather than a
# number someone chose.
#
# WHY THE CONTRACT STILL DECLARES THEM RATHER THAN ONLY EMITTING THEM. The
# consumer is a Fortran namelist. `landmod_nl` takes scalars; a run RECORDS the
# pair it integrated, which is what makes the soil heat solver's endpoints an
# arm rather than a compiled accident; and `landmod.f90` carries the same pair as
# a compiled default so a run that declares nothing reproduces it. A namelist
# writer that reached into this report instead would make a run's arm depend on
# whichever soil report happened to be on disk, unnamespaced and unstamped. What
# changes here is not that the scalar stops being declared. It is that nobody
# types it.
#
# THE THERMAL ENDPOINTS ARE NOT IN THIS SET and are a different class. They are
# evaluated at a DECLARED column -- `endpoints.evaluated_at`'s bulk density and
# quartz fraction -- so they move only when that column is re-declared, and the
# check that the declared column is still the soil map's median is separate and
# is the one that fires when it stops being.
DERIVED_MAPPING_SCALARS = ("sr_at_wilting_point", "sr_at_field_capacity")
DERIVED_MAPPING_PAIRS = ("sr_at_wilting_point_p5_p95",
                         "sr_at_field_capacity_p5_p95")


def _rewrite_key(text: str, key: str, rendered: str) -> tuple[str, str | None]:
    """Replace one `key: value` line in the declaration, keeping its indent.

    Line-targeted rather than a YAML round trip because the declaration is
    mostly argument: a dump would lose every comment in it, and the comments are
    the part a reader needs.
    """
    pattern = re.compile(rf"^(?P<indent>[ \t]*){re.escape(key)}:[ \t]*"
                         r"(?P<old>[^\n#]*?)[ \t]*$", re.MULTILINE)
    found = pattern.findall(text)
    if len(found) != 1:
        raise SystemExit(
            f"{DECLARATION.name} states `{key}:` {len(found)} times and this "
            "rewrite needs exactly one. The declaration has been restructured "
            "and `--update-declaration` has to be restructured with it.")
    match = pattern.search(text)
    old = match.group("old")
    if old == rendered:
        return text, None
    return (text[:match.start()] + match.group("indent") + key + ": " + rendered
            + text[match.end():], f"{key}: {old} -> {rendered}")


def update_declaration(measured: dict) -> list[str]:
    """Write the derived saturation mapping into the declaration.

    Returns one line per number that moved, empty when the declaration already
    states the measurement. Writes the file only when something moved.
    """
    text = DECLARATION.read_text(encoding="utf-8")
    changed: list[str] = []
    for key in DERIVED_MAPPING_SCALARS:
        text, note = _rewrite_key(text, key, f"{measured[key]:.4f}")
        if note:
            changed.append(note)
    for key in DERIVED_MAPPING_PAIRS:
        low, high = measured[key]
        text, note = _rewrite_key(text, key, f"[{low:.4f}, {high:.4f}]")
        if note:
            changed.append(note)
    if changed:
        DECLARATION.write_text(text, encoding="utf-8")
    return changed


def _fortran_default(name: str) -> float | None:
    """One `real :: name = value` default out of `landmod.f90`, or None."""
    if not LANDMOD_SOURCE.is_file():
        return None
    text = LANDMOD_SOURCE.read_text(encoding="utf-8", errors="replace")
    m = re.search(rf"^[ \t]*real[ \t]*::[ \t]*{name}[ \t]*=[ \t]*"
                  r"([-+0-9.eEdD]+)", text, re.IGNORECASE | re.MULTILINE)
    if not m:
        return None
    return float(m.group(1).replace("d", "e").replace("D", "e"))


def scan_restatements(decl: dict, config: dict) -> list[dict]:
    """Every other place the declared saturation endpoints are written.

    The declaration is one number and a run reaches it through copies: two
    `landmod_nl` compiled defaults for the soil heat solver, one more for the
    surface layer the moisture-dependent albedo mixes on, and the project config
    keys `run_exoplasim.py` writes into the namelist. Each row says what that
    copy states and what the contract says it should, so a rebuild that moves the
    median produces the work list instead of leaving it to be found.

    Compared against the DECLARATION and not against the measurement, because
    that is the fact they restate; that the declaration IS the measurement is
    `check_thermal_against_states`'s job and is checked separately.
    """
    mapping = decl["thermal"]["saturation_mapping"]
    tol = float(mapping["tolerance"])
    wp = float(mapping["sr_at_wilting_point"])
    fc = float(mapping["sr_at_field_capacity"])
    surface = decl["surface_layer"]["saturation_mapping"]
    node = decl
    for part in str(surface["sr_at_field_capacity_from"]).split("."):
        node = node[part]
    skin_fc = float(node)

    surf_cfg = config.get("surface", {}) or {}
    thermal_cfg = surf_cfg.get("soil_thermal", {}) or {}
    albedo_cfg = surf_cfg.get("soil_albedo_moisture", {}) or {}

    def stated(block: dict, key: str) -> float | None:
        value = block.get(key)
        return None if value is None else float(value)

    landmod = "vendor/exoplasim/exoplasim/plasim/src/landmod.f90"
    by_import = "exoplasim/scripts/run_exoplasim.py, at import"
    by_gate = "scripts/check_consistency.py"
    rows = [
        {"where": landmod, "key": "soilsrwp", "declared": wp,
         "states": _fortran_default("soilsrwp"), "enforced_by": by_import},
        {"where": landmod, "key": "soilsrfc", "declared": fc,
         "states": _fortran_default("soilsrfc"), "enforced_by": by_import},
        {"where": landmod, "key": "skinsrfc", "declared": skin_fc,
         "states": _fortran_default("skinsrfc"), "enforced_by": by_import},
        {"where": "config/planet.yaml",
         "key": "surface.soil_albedo_moisture.saturation_at_full_layer",
         "declared": skin_fc,
         "states": stated(albedo_cfg, "saturation_at_full_layer"),
         "enforced_by": by_gate},
        {"where": "config/planet.yaml",
         "key": "surface.soil_thermal.saturation_at_empty_store",
         "declared": wp,
         "states": stated(thermal_cfg, "saturation_at_empty_store"),
         "enforced_by": by_gate},
        {"where": "config/planet.yaml",
         "key": "surface.soil_thermal.saturation_at_full_store",
         "declared": fc,
         "states": stated(thermal_cfg, "saturation_at_full_store"),
         "enforced_by": by_gate},
    ]
    for row in rows:
        if row["states"] is None:
            # UNSET IS A LEGITIMATE STATE AND NOT A GAP. `run_exoplasim.py`
            # writes every one of these keys into `landmod_nl` whatever the
            # config says; a key the config leaves out is written at
            # `landmod.f90`'s compiled default, which is itself a row above and
            # is checked against this contract. So an unset key is one copy of
            # the number instead of two, and adding the key would add a
            # restatement rather than remove one.
            row["status"] = "unset"
        elif abs(row["states"] - row["declared"]) > tol:
            row["status"] = "stale"
        else:
            row["status"] = "agrees"
    return rows


# ---------------------------------------------------------------------------
# 1. The declaration
# ---------------------------------------------------------------------------

def check_thermal_against_states(decl: dict, soil: dict[str, np.ndarray],
                                 gravity_m_s2: float) -> tuple[dict, list[str]]:
    """The declared thermal reduction, recomputed from this build's own states.

    `thermal.saturation_mapping` and `thermal.endpoints` are SCALARS standing in
    for per-cell quantities this contract already emits, and a scalar standing in
    for a distribution is honest only while it is the statistic it says it is.
    So both are recomputed here and the contract is refused when either has
    drifted: the mapping's two endpoints against the medians of theta_wp/theta_s
    and theta_fc/theta_s, and the conductivity and heat capacity endpoints
    against Johansen's relations at the declared column.

    THE RELATIONS COME FROM `analysis/soil_thermal_inertia.py` and are not
    restated here. That module is where Farouki's coefficients and his two
    texture branches live, and a second copy would be free to agree with
    neither.
    """
    sys.path.insert(0, str(PROJECT_ROOT / "analysis"))
    import soil_thermal_inertia as sti

    _, theta_s, theta_fc, theta_wp = cosby_states_from_contract(
        soil, decl, gravity_m_s2)
    sr_wp = theta_wp / theta_s
    sr_fc = theta_fc / theta_s

    mapping = decl["thermal"]["saturation_mapping"]
    ends = decl["thermal"]["endpoints"]
    rho_b = float(ends["evaluated_at"]["bulk_density_kg_m3"])
    quartz = float(ends["evaluated_at"]["quartz_fraction"])

    measured = {
        "sr_at_wilting_point": float(np.median(sr_wp)),
        "sr_at_field_capacity": float(np.median(sr_fc)),
        "sr_at_wilting_point_p5_p95": [float(np.percentile(sr_wp, 5)),
                                       float(np.percentile(sr_wp, 95))],
        "sr_at_field_capacity_p5_p95": [float(np.percentile(sr_fc, 5)),
                                        float(np.percentile(sr_fc, 95))],
        "median_bulk_density_kg_m3": float(np.median(soil["bulkdensity"])),
        "dry_conductivity_w_m_k": float(sti.conductivity(0.0, quartz, rho_b)),
        "saturated_conductivity_w_m_k": float(sti.conductivity(1.0, quartz, rho_b)),
        "dry_heat_capacity_j_m3_k": float(sti.heat_capacity(0.0, rho_b)),
        "saturated_heat_capacity_j_m3_k": float(sti.heat_capacity(1.0, rho_b)),
    }

    # What the store can actually reach, which is the point of the whole block:
    # an empty store is the wilting point, not a dry soil.
    reach = {}
    for label, sr in (("empty_store", measured["sr_at_wilting_point"]),
                      ("full_store", measured["sr_at_field_capacity"]),
                      ("air_dry", sti.SATURATION_RANGE[0]),
                      ("saturated", 1.0)):
        reach[label] = round(float(sti.inertia(sr, quartz, rho_b)), 1)
    reach["store_range_factor"] = round(reach["full_store"] / reach["empty_store"], 3)
    reach["full_sweep_factor"] = round(reach["saturated"] / reach["air_dry"], 3)
    measured["thermal_inertia_j_m2_k_s05"] = reach

    bad: list[str] = []
    tol_map = float(mapping["tolerance"])
    for key in ("sr_at_wilting_point", "sr_at_field_capacity"):
        if abs(float(mapping[key]) - measured[key]) > tol_map:
            bad.append(
                f"thermal.saturation_mapping.{key} declares {mapping[key]} and "
                f"this build's states give {measured[key]:.4f}. It is declared "
                f"as the median of a distribution this contract emits, so a "
                f"drift past {tol_map} means it is no longer that statistic")
    tol_end = float(ends["tolerance"])
    for path, declared_value, name in (
            (("dry", "thermal_conductivity_w_m_k"),
             ends["dry"]["thermal_conductivity_w_m_k"], "dry_conductivity_w_m_k"),
            (("saturated", "thermal_conductivity_w_m_k"),
             ends["saturated"]["thermal_conductivity_w_m_k"],
             "saturated_conductivity_w_m_k"),
            (("dry", "volumetric_heat_capacity_j_m3_k"),
             ends["dry"]["volumetric_heat_capacity_j_m3_k"],
             "dry_heat_capacity_j_m3_k"),
            (("saturated", "volumetric_heat_capacity_j_m3_k"),
             ends["saturated"]["volumetric_heat_capacity_j_m3_k"],
             "saturated_heat_capacity_j_m3_k")):
        got = measured[name]
        if abs(float(declared_value) - got) > tol_end * abs(got):
            bad.append(
                f"thermal.endpoints.{'.'.join(path)} declares {declared_value} "
                f"and Johansen's relations at the declared column give {got:.6g}. "
                f"The endpoint is derived, not chosen, so the two must agree")
    if abs(measured["median_bulk_density_kg_m3"] - rho_b) > 1.0:
        bad.append(
            f"thermal.endpoints.evaluated_at.bulk_density_kg_m3 is {rho_b} and "
            f"the soil map's median is "
            f"{measured['median_bulk_density_kg_m3']:.1f}. The endpoints are "
            "evaluated at the build's own median column and that is the "
            "statistic they are declared as")
    return measured, bad


def cosby_states_from_contract(soil: dict[str, np.ndarray], decl: dict,
                               gravity_m_s2: float):
    """The adopted per-cell states, in the one place they are derived."""
    return contract_states(soil["sand"], soil["clay"], decl, gravity_m_s2)


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

    for section in ("states", "flow", "thermal", "materials", "uncertainty",
                    "surface_layer"):
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
    """Pedology's declared endmember capacity against the shipped Cosby inversion.

    THE TWO DERIVATIONS THIS REPLACED, kept because the attribution is what
    made the replacement legible. `soilmap.txt`'s `awc` is pedology's endmember
    mixture with its additive organic and allophane terms; `lpj_capacity_mm` is
    the Cosby retention inversion at Cosby's own suction with the regolith and
    bedrock rescaling, transcribed from the source as it stood. Neither is what
    a consumer installs now: both read the emitted states.
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
        "exoplasim_source": "what dwmax WAS: the pedology awc column converted "
                            "mm to m. It is now the emitted states file's "
                            "awc_mm column on the same terms",
        "lpj_over_pedology_ratio": percentiles(ratio),
        "absolute_difference_mm": percentiles(difference),
        "pearson_correlation": correlation,
        "reading": "The two fields correlate because both start from related "
                   "texture and depth. They are not the same capacity: one is a "
                   "declared endmember mixture with additive organic and "
                   "allophane terms, the other is a Cosby retention inversion "
                   "with no andic term and no bulk density. That is the "
                   "disagreement the adopted case ended.",
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

    `contract_states` clamps both defining pressures at the air-entry pressure,
    so a cell inside that region does not return a field capacity above
    saturation; it returns saturation, which is the physically right answer and
    is a DIFFERENT number from the one the retention curve would give. The map
    is still measured against the line as a MARGIN, because a count of zero
    says nothing about how close the nearest cell is, and because a map that
    moved inside would have its field capacity set by the clamp rather than by
    its texture.

    The line is derived at Cosby's own field-capacity suction, which is the
    boundary the regression coefficients give in closed form. This world's
    field-capacity pressure is LARGER, so the region where the clamp binds is
    strictly smaller than this line and the margin reported here is the
    conservative one.
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
                   "WORLD-NGA10 clamped both heads at the air-entry value in "
                   "get_mineral; WORLD-OF6N moved the derivation into "
                   "contract_states and the clamp with it, so the closure "
                   "cannot return a field capacity above saturation for any "
                   "texture. No cell of this map is inside the region, so the "
                   "clamp changes no current soil property; the margin is what "
                   "says how far that is from being true by accident.",
        "bound_is_conservative": "the line is derived at Cosby's own "
                                 "field-capacity suction. The adopted suction "
                                 "is larger, so the region where the clamp "
                                 "binds is strictly smaller than this",
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
        "applied_by_any_consumer": True,
        "applied_note": "APPLIED. WORLD-OF6N adopted the contract read literally, "
                        "so the emitted states evaluate field capacity at "
                        "rho_w * g * L with this planet's gravity and every "
                        "consumer installs them. Pedology's declared endmember "
                        "capacities are Earth field observations and are no "
                        "longer what any consumer reads.",
        "bracket_note": "an UPPER BOUND on the size of the effect and adopted as "
                        "one: it holds the Clapp-Hornberger exponent fixed while "
                        "gravity moves, so it is the whole of the shift a "
                        "gravity-invariant curve can produce, and the "
                        "within-texture-class spread Cosby reports is larger "
                        "than the shift -- `uncertainty_cases` measures both "
                        "moves in the same capacity. A declared limitation of the adopted "
                        "value, not a reason to leave the states at Earth's "
                        "suction: the shift moves every cell the same way, so "
                        "unlike a scatter it does not average out over the map.",
    }


# ---------------------------------------------------------------------------
# 5. The adopted central case, and what it moved
# ---------------------------------------------------------------------------

def adopted_case(decl: dict, soil: dict[str, np.ndarray],
                 gravity_m_s2: float) -> dict:
    """The one derivation both consumers read, and what adopting it moved.

    WORLD-OF6N is settled. The central case is the contract read literally: the
    named closure -- Clapp-Hornberger on Cosby Table 4 -- evaluated at the
    suction the contract's own potential convention declares, `rho_w * g * L`
    with this planet's gravity. WORLD-SLPA is the same change reached from the
    gravity side and is subsumed by it, because landing both would apply the
    gravity conversion once and a frame-mixing artifact once.

    The two alternatives are kept as what was NOT adopted and what each would
    have cost, which is the attribution the removal needed. Neither is an open
    option.

    THE COMPARISON'S BASELINE MOVED WITH WORLD-7702. Pedology's `awc` column
    now stops at the declared column base, so the ratios below are against a
    capacity that no longer follows regolith past it. That was the larger of
    the two disagreements between the consumers, and the adopted case settles
    it by construction: the adopted capacity is an integral over the declared
    layers, so both consumers cover the same column.
    """
    pedology_mm = soil["awc"]
    lpj_shipped_mm = lpj_capacity_mm(soil["sand"], soil["clay"],
                                     soil["depth"], soil["bedrockfrac"])
    _, _, theta_fc, theta_wp = contract_states(soil["sand"], soil["clay"],
                                               decl, gravity_m_s2)
    usable = layer_usable_fraction(soil["depth"], soil["bedrockfrac"], decl)
    adopted_mm = column_capacity_mm(theta_fc, theta_wp, usable, decl)

    def moves(new, old):
        ratio = new / np.maximum(old, 1e-12)
        return {"median_mm_before": float(np.median(old)),
                "median_mm_after": float(np.median(new)),
                "land_mean_mm_before": float(old.mean()),
                "land_mean_mm_after": float(new.mean()),
                "ratio": percentiles(ratio)}

    adopted = {
        "name": "closure_at_vesper_suction",
        "what": "the contract's named closure at the suction its own potential "
                "convention declares, rho_w * g * L at this planet's gravity. "
                "The closure is named and field capacity is declared as a "
                "drainage equilibrium over a declared length, so this is the "
                "contract read literally and no part of it is a new choice",
        "central_case_mm": percentiles(adopted_mm),
        "exoplasim_dwmax": moves(adopted_mm, pedology_mm),
        "lpj_guess_capacity": moves(adopted_mm, lpj_shipped_mm),
        "deleted": "pedology's endmember volumetric capacities, its additive "
                   "organic term and its allophane term as the source of any "
                   "consumer's capacity; soilinput.cpp's get_mineral Cosby "
                   "inversion; and vesperinput.cpp's apply_regolith_depth "
                   "rescaling, whose rule is now layer_usable_fraction here",
        "limitation": "the gravity shift holds the Clapp-Hornberger exponent "
                      "fixed while gravity moves, so the adopted field capacity "
                      "is an UPPER BOUND on how far a gravity-invariant "
                      "retention curve can fall.",
        "uncertainty": "Cosby's within-texture-class residual spread, carried in "
                       "`uncertainty_cases` and declared in the contract's "
                       "uncertainty block. It is larger than the gravity shift, "
                       "and that comparison is measured there rather than "
                       "asserted here.",
    }

    not_adopted = {
        "closure_at_earth_suction": {
            "what": "the same closure evaluated at Cosby's own field-capacity "
                    "suction, which is the pressure the shipped LPJ-GUESS code "
                    "used. Not adopted: it reads the contract's closure and "
                    "ignores the contract's potential convention, which "
                    "declares field capacity as a drainage equilibrium and "
                    "therefore at this world's gravity",
            "central_case_mm": percentiles(lpj_shipped_mm),
            "would_have_moved_exoplasim": moves(lpj_shipped_mm, pedology_mm),
            "would_have_moved_lpj_guess": moves(lpj_shipped_mm, lpj_shipped_mm),
        },
        "pedology_endmember_mixture": {
            "what": "pedology's declared endmember volumetric capacities with "
                    "its additive organic and allophane terms. Not adopted, "
                    "and not a symmetric option: it is not a retention curve "
                    "at all -- no saturation, no air entry, no exponent -- so "
                    "taking it would have replaced the contract's named "
                    "closure and left saturation and matric potential with "
                    "nowhere to come from",
            "central_case_mm": percentiles(pedology_mm),
            "would_have_moved_exoplasim": moves(pedology_mm, pedology_mm),
            "would_have_moved_lpj_guess": moves(pedology_mm, lpj_shipped_mm),
        },
    }

    andic_mm = (float(_pedogenesis()["andisol"]["volumetric_capacity_allophane"])
                * soil["andic"]
                * np.minimum(soil["depth"],
                             float(decl["geometry"]["column_base_m"]))
                * 1000.0)
    return {
        "adopted": adopted,
        "not_adopted": not_adopted,
        "andic_term_mm": percentiles(andic_mm),
        "andic_cells": int((soil["andic"] > 0.0).sum()),
        "andic_share_of_pedology_awc": percentiles(
            andic_mm / np.maximum(pedology_mm, 1e-12)),
        "andic_is_not_the_disagreement": (
            "the allophane term LPJ-GUESS has no equivalent of contributes a "
            "median of nothing and a p90 of a few per cent of pedology's "
            "capacity, against a median absolute gap between the two "
            "derivations of tens of millimetres. What separated them is that "
            "the declared endmember capacities and the Cosby retention "
            "inversion are different numbers for the same texture, not that "
            "one carries a material the other does not. The adopted closure "
            "carries no allophane term either, so the andic hydraulic effect "
            "is now a DECLARED absence rather than a disagreement between two "
            "consumers"),
        "median_absolute_gap_mm": float(np.median(
            np.abs(lpj_shipped_mm - pedology_mm))),
        "not_measurable_here": (
            "the change in land runoff ratio and in E over R. Both need a "
            "paired baseline, pedology/scripts/probe_runoff_response.py needs "
            "a climatology to back-solve potential evaporation from, and "
            "config/planet.yaml declares baseline_climatology null. On this "
            "build the capacity change is costed and its climate consequence "
            "is not available at all"),
        "not_execution_verified": (
            "nothing on the LPJ-GUESS side. It does not build on this tree "
            "until a baseline climatology exists, so every LPJ number here is "
            "a transcription of the vendored source evaluated offline"),
        "owner": "world-of6n",
    }


# ---------------------------------------------------------------------------
# 6. The uncertainty of the adopted case
# ---------------------------------------------------------------------------

# The corners of the one-sigma box in Cosby's three parameters. Three
# parameters and two signs is eight cases, all of them evaluated, because
# Cosby reports no correlation between the three residuals and so nothing here
# can say which corner a real soil sits nearest. Taking the extremes over the
# corners is therefore an ENVELOPE and is labelled as one everywhere it is
# reported: it is what the spread can be worth, not a percentile of a
# distribution.
CASE_DEVIATE = 1.0
CASE_CORNERS = tuple((zb, zp, zt)
                     for zb in (-CASE_DEVIATE, CASE_DEVIATE)
                     for zp in (-CASE_DEVIATE, CASE_DEVIATE)
                     for zt in (-CASE_DEVIATE, CASE_DEVIATE))


def uncertainty_cases(decl: dict, soil: dict[str, np.ndarray],
                      gravity_m_s2: float) -> dict:
    """What the adopted states are worth, from Cosby's within-class spread.

    THE CENTRAL CASE'S OWN UNCERTAINTY. Texture explains the mean structure of
    the retention parameters over 1,448 samples, and the parameter spread that
    survives inside a texture class is large -- Cosby's Table 3 puts the
    exponent's within-class standard deviation at a third of its mean. Every
    consumer installs the central case, so without this the states carry no
    declared spread at all and the percentiles in this report are a spread
    ACROSS THE MAP, which is a different fact about a different thing.

    Applied COHERENTLY, which is what makes this a bracket rather than a
    smoothing: one set of deviates moves the parameters, the parameters set
    every state, and the states apply to all fifteen layers of the cell. Drawing
    a layer at a time would average the spread away and report a precision the
    contract does not have.

    The low and high cases are the per-cell extremes over the corners of the
    one-sigma box, which is a bound and not a quantile: Cosby reports the three
    residual spreads separately and no correlation between them, so a case that
    takes a low exponent with a high porosity cannot be ruled out from the
    published table, and nothing here can weight it either.
    """
    sand, clay = soil["sand"], soil["clay"]
    usable = layer_usable_fraction(soil["depth"], soil["bedrockfrac"], decl)

    def case(deviates):
        b, theta_s, theta_fc, theta_wp = contract_states(
            sand, clay, decl, gravity_m_s2, deviates)
        return {"b": b, "theta_s": theta_s, "theta_fc": theta_fc,
                "theta_wp": theta_wp,
                "awc_mm": column_capacity_mm(theta_fc, theta_wp, usable, decl)}

    central = case(None)
    corners = [case(z) for z in CASE_CORNERS]
    names = ("b", "theta_s", "theta_fc", "theta_wp", "awc_mm")
    low = {name: np.min([c[name] for c in corners], axis=0) for name in names}
    high = {name: np.max([c[name] for c in corners], axis=0) for name in names}

    # The ordering the contract declares has to survive every case, or a case
    # is not a soil. This is the check with a right answer.
    violations = 0
    for c in corners:
        violations += int(np.count_nonzero(
            (c["theta_wp"] > c["theta_fc"] + 1e-12)
            | (c["theta_fc"] > c["theta_s"] + 1e-12)))

    driving = np.argmax([c["awc_mm"] for c in corners], axis=0)
    spread = cosby_parameter_spread(sand, clay)

    # The contract declares the gravity shift in field capacity as an upper
    # bound and says this spread is larger than it. Both are moves in the same
    # capacity, so the comparison is arithmetic rather than assertion: the
    # gravity move is the central case against the same closure at Earth's
    # gravity, and the case move is the wider side of the envelope.
    earth_g = float(decl["potential_convention"]["earth_gravity_m_s2"])
    _, _, fc_earth, wp_earth = contract_states(sand, clay, decl, earth_g)
    awc_earth = column_capacity_mm(fc_earth, wp_earth, usable, decl)
    gravity_move = np.abs(central["awc_mm"] / np.maximum(awc_earth, 1e-12) - 1.0)
    case_move = np.maximum(
        np.abs(low["awc_mm"] / central["awc_mm"] - 1.0),
        np.abs(high["awc_mm"] / central["awc_mm"] - 1.0))
    against_gravity = {
        "gravity_move_in_awc": percentiles(gravity_move),
        "case_move_in_awc": percentiles(case_move),
        "case_over_gravity": percentiles(case_move / np.maximum(gravity_move,
                                                                1e-12)),
        "cells_where_the_case_move_is_smaller": int(
            np.count_nonzero(case_move < gravity_move)),
        "note": "the gravity shift is one-signed and moves every cell the same "
                "way, so it does not average out over the map; this spread is a "
                "scatter and does. They are not interchangeable and the "
                "comparison is of magnitude only.",
    }
    return {
        "source": "cosby1984 Table 4, the S.D. rows: the within-texture-class "
                  "residual spread of the same parameters whose means are the "
                  "adopted central case",
        "reference": "references/cosby_1984_a-statistical-exploration-of-the-"
                     "relationships-of-soil-moisture-charac.pdf",
        "deviate": CASE_DEVIATE,
        "corners": len(CASE_CORNERS),
        "applied": decl["uncertainty"]["applied"],
        "is_an_envelope": (
            "low and high are the per-cell extremes over the corners of the "
            "one-sigma box in b, log Psi_s and Theta_s. Cosby reports the three "
            "residual spreads separately and no correlation between them, so "
            "this is a bound on what the spread can be worth and not a "
            "percentile of a distribution. A correlation, if it were published, "
            "could only narrow it."),
        "within_class_sd": {
            "b": percentiles(spread["b"]),
            "log_psi_s": percentiles(spread["log_psi_s"]),
            "theta_s_volumetric": percentiles(spread["theta_s"]),
            "b_relative_to_central": percentiles(spread["b"] / central["b"]),
            "theta_s_relative_to_central": percentiles(
                spread["theta_s"] / central["theta_s"]),
        },
        "awc_mm": {
            "central": percentiles(central["awc_mm"]),
            "low": percentiles(low["awc_mm"]),
            "high": percentiles(high["awc_mm"]),
            "low_over_central": percentiles(low["awc_mm"] / central["awc_mm"]),
            "high_over_central": percentiles(high["awc_mm"] / central["awc_mm"]),
        },
        "states": {
            name: {"central": percentiles(central[name]),
                   "low": percentiles(low[name]),
                   "high": percentiles(high[name])}
            for name in ("theta_s", "theta_fc", "theta_wp", "b")
        },
        "ordering_violations_over_all_cases": violations,
        "cells_by_high_corner": {
            str(CASE_CORNERS[corner]): int(np.count_nonzero(driving == corner))
            for corner in range(len(CASE_CORNERS))
            if int(np.count_nonzero(driving == corner))
        },
        "installed_by_any_consumer": False,
        "installed_note": "the emitted states are the central case and every "
                          "consumer installs that. A case is a states file "
                          "away: `contract_states` takes the deviates and "
                          "nothing downstream of it knows which case it read.",
        "against_the_gravity_shift": against_gravity,
    }


# ---------------------------------------------------------------------------
# 7. The emission. What the consumers read instead of deriving.
# ---------------------------------------------------------------------------

# Shared with pedology/scripts/build_soil.py and
# biosphere/scripts/build_lpj_driver.py: every one of these files is keyed on an
# exactly-compared pair of doubles, so they must all round identically or a
# lookup silently misses.
COORD_DECIMALS = 4

STATE_COLUMNS = ("b", "theta_s", "theta_fc", "theta_wp", "awc_mm",
                 "evaporable_mm")


def emit_states(decl: dict, soil: dict[str, np.ndarray], gravity_m_s2: float,
                path: Path) -> dict:
    """Write the contract's per-cell states, which is what a consumer reads.

    One row per land cell of the soil map, in the soil map's own order and on
    its own coordinates. The columns are the closure's exponent, the three
    volumetric states, the plant-available capacity of the whole physical
    column, and the per-layer usable fraction the weathered-bedrock rule gives.

    WHY THE PER-LAYER FRACTION AND NOT PER-LAYER CAPACITIES. The same fraction
    scales a layer's available water, its wilting point and its saturation,
    because material below the regolith contact holds less of everything. One
    number per layer instead of three is what stops the three drifting apart in
    a consumer, and it is what makes the column capacity here and the sum over
    a consumer's layers the same arithmetic rather than two that agree today.

    TWO CAPACITY COLUMNS, AND THEY ARE TWO QUANTITIES. `awc_mm` is
    plant-available over the whole physical column and is what LPJ-GUESS
    installs; `evaporable_mm` adds the surface layer's water between air dry
    and the wilting point and is what ExoPlaSim installs as `dwmax`. A root
    cannot reach that water and a bare drying surface can, so this is not one
    number written twice. `surface_layer_report` carries the difference and
    checks the identity that defines it. `soilmap.txt`'s `awc` is a third
    thing: pedology's own declared endmember mixture, which this report costs
    against the adopted closure and which no consumer reads.
    """
    b, theta_s, theta_fc, theta_wp = contract_states(
        soil["sand"], soil["clay"], decl, gravity_m_s2)
    usable = layer_usable_fraction(soil["depth"], soil["bedrockfrac"], decl)
    capacity = column_capacity_mm(theta_fc, theta_wp, usable, decl)
    evaporable = evaporable_capacity_mm(theta_fc, theta_wp, usable, decl)
    nlayer = int(decl["geometry"]["physical_layer_count"])

    header = ("Lon Lat " + " ".join(STATE_COLUMNS) + " "
              + " ".join(f"u{layer:02d}" for layer in range(nlayer)))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(header + "\n")
        for cell in range(capacity.size):
            row = [f"{soil['Lon'][cell]:.{COORD_DECIMALS}f}",
                   f"{soil['Lat'][cell]:.{COORD_DECIMALS}f}",
                   f"{b[cell]:.5f}", f"{theta_s[cell]:.6f}",
                   f"{theta_fc[cell]:.6f}", f"{theta_wp[cell]:.6f}",
                   f"{capacity[cell]:.3f}", f"{evaporable[cell]:.3f}"]
            row.extend(f"{usable[cell, layer]:.6f}" for layer in range(nlayer))
            handle.write(" ".join(row) + "\n")

    return {
        "path": str(path.resolve().relative_to(PROJECT_ROOT)),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "cells": int(capacity.size),
        "physical_layers": nlayer,
        "gravity_m_s2": gravity_m_s2,
        "awc_mm": percentiles(capacity),
        "evaporable_mm": percentiles(evaporable),
        "saturation_volumetric": percentiles(theta_s),
        "field_capacity_volumetric": percentiles(theta_fc),
        "wilting_point_volumetric": percentiles(theta_wp),
        "closure": decl["retention_closure"]["name"],
        "what_it_is": "the adopted per-cell retention states and the "
                      "weathered-bedrock usable fraction per physical layer. "
                      "ExoPlaSim installs evaporable_mm as dwmax, LPJ-GUESS "
                      "takes awc_mm and the states and the fractions through "
                      "the driver file; neither derives them.",
    }


def read_states(path: Path) -> dict[str, np.ndarray]:
    """Read an emitted states file. The one reader, shared by every consumer."""
    data = np.genfromtxt(path, names=True)
    return {name: np.atleast_1d(np.asarray(data[name], dtype=float))
            for name in data.dtype.names}


def usable_columns(states: dict[str, np.ndarray]) -> np.ndarray:
    """The per-layer usable fractions of an emitted states file, as one array."""
    names = sorted(name for name in states if name.startswith("u")
                   and name[1:].isdigit())
    return np.stack([states[name] for name in names], axis=1)


def _pedogenesis() -> dict:
    return yaml.safe_load((COMPONENT_ROOT / "config" / "pedogenesis.yaml")
                          .read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------

def surface_layer_report(decl: dict, soil: dict[str, np.ndarray],
                        gravity_m_s2: float) -> tuple[dict, list[str]]:
    """What the surface layer costs the column, and the identity that defines it.

    Three things, and all three are numbers rather than claims: how much water
    the column gained by cutting its top 0.02 m from air dry instead of from
    the wilting point, how large that is against the capacity it is added to,
    and the check that the two emitted capacity columns differ by exactly that
    and by nothing else.
    """
    _b, _theta_s, theta_fc, theta_wp = contract_states(
        soil["sand"], soil["clay"], decl, gravity_m_s2)
    usable = layer_usable_fraction(soil["depth"], soil["bedrockfrac"], decl)
    awc = column_capacity_mm(theta_fc, theta_wp, usable, decl)
    cap, increment = surface_layer_capacity_mm(theta_fc, theta_wp, usable, decl)
    evaporable = evaporable_capacity_mm(theta_fc, theta_wp, usable, decl)

    residual = float(np.max(np.abs((evaporable - awc) - increment)))
    bad: list[str] = []
    if residual > 1.0e-9:
        bad.append(
            f"evaporable_mm - awc_mm differs from the surface layer's "
            f"sub-wilting increment by up to {residual:.3e} mm. The two "
            "capacity columns are supposed to differ by exactly that water and "
            "by nothing else, so one of them is being built from something the "
            "other does not carry.")
    surf = decl["surface_layer"]
    split = check_split_invariance()
    if not split["passes"]:
        bad.append(
            "splitting the top water layer moved water: the two- and "
            "three-layer columns disagree by "
            f"{split['worst_total_disagreement_mm']:.3e} mm on the total or "
            f"{split['worst_runoff_disagreement_mm']:.3e} mm on the runoff. "
            "The surface layer is only free because that difference is zero.")
    return {
        "thickness_m": float(surf["thickness_m"]),
        "air_dry_water_content": float(surf["air_dry_water_content"]),
        "capacity_mm": percentiles(cap),
        "sub_wilting_increment_mm": percentiles(increment),
        "awc_mm": percentiles(awc),
        "evaporable_mm": percentiles(evaporable),
        "increment_over_awc": percentiles(increment / np.maximum(awc, 1.0e-12)),
        "air_dry_overdraw_mm": {
            "what": "the water the surface layer can give up that a real soil "
                    "would hold by adsorption at the ambient humidity. "
                    "`air_dry_water_content` is declared zero, so the overdraw "
                    "is the sub-wilting increment times theta_ad over theta_wp, "
                    "an unknown fraction under one; the BOUND and its direction "
                    "are what is known",
            "bounded_by_mm": percentiles(increment),
        },
        "identity_residual_mm": residual,
        "split_invariance": split,
        "saturation_endpoints": {
            "sr_at_air_dry": float(surf["saturation_mapping"]["sr_at_air_dry"]),
            "sr_at_field_capacity": float(
                decl["thermal"]["saturation_mapping"]["sr_at_field_capacity"]),
            "note": "the surface layer's own mapping. It differs from the "
                    "thermal one in the lower endpoint only, because its "
                    "capacity is cut from air dry and the layers below it are "
                    "cut from the wilting point.",
        },
        "consumers": surf["consumers"],
    }, bad


def build_report(decl: dict, config: dict, soil_map: Path,
                 gravity_m_s2: float,
                 states_path: Path | None = None) -> dict:
    soil = read_soil_map(soil_map)
    failures = check_declaration(decl)
    mutations, defects = run_declaration_mutations()
    frames = frame_consistency(soil, decl, gravity_m_s2)
    defects = defects + frames["defects"]
    thermal_measured, thermal_bad = check_thermal_against_states(
        decl, soil, gravity_m_s2)
    failures = failures + thermal_bad
    surface, surface_bad = surface_layer_report(decl, soil, gravity_m_s2)
    failures = failures + surface_bad
    emitted = (emit_states(decl, soil, gravity_m_s2, states_path)
               if states_path is not None else None)
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
        "frame_consistency": frames,
        "emitted_states": emitted,
        "adopted_case": adopted_case(decl, soil, gravity_m_s2),
        "uncertainty_cases": uncertainty_cases(decl, soil, gravity_m_s2),
        "undeclared_properties": undeclared_properties(decl),
        "thermal_reduction": thermal_measured,
        "surface_layer": surface,
        "restatements": scan_restatements(decl, config),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--soil-map", type=Path, default=None,
                        help="defaults to the configured build's soil map")
    parser.add_argument("--strict", action="store_true",
                        help="exit 1 while any contract property is still undeclared")
    parser.add_argument("--states", type=Path, default=None,
                        help="where to write the per-cell states every consumer "
                             "reads; defaults to the configured build's")
    parser.add_argument("--no-emit", action="store_true",
                        help="check and report only. The states file is what "
                             "ExoPlaSim and LPJ-GUESS install, so not writing "
                             "it leaves whatever is on disk in place")
    parser.add_argument("--update-declaration", action="store_true",
                        help="write the derived saturation mapping into "
                             "pedology/config/land_column_properties.yaml from "
                             "this build's own states, instead of refusing "
                             "because a hand-typed copy of it has drifted")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    gravity = float(config["planet"]["gravity_m_s2"])
    sys.path.insert(0, str(PROJECT_ROOT / "lib"))
    import builds
    if args.soil_map is None:
        args.soil_map = builds.soilmap(config)
    if args.states is None:
        args.states = builds.land_column_states(config)
    # THE TWO PATHS ARE PER RUNG AND NOTHING TIED THEM TOGETHER. This script is
    # a pure carrier -- the states come from the soil map's rows and no climate
    # -- so it already runs at any rung once a soil map exists there, and both
    # paths are overridable. That is exactly what makes the mismatch reachable:
    # `--soil-map soilmap_T42.txt` with the default `--states` writes T42 rows
    # into `land_column_states_T21.txt`, and the consumers key the file by rung.
    # `build_surface_soil_water.py` would then read a states file whose cell
    # count is not its grid's. Both names carry their rung, so this is an
    # identity check and not a guess; a path that does not name a rung is left
    # alone, because a caller writing to a scratch file means it.
    import rungs
    def _named_rung(path: Path) -> str | None:
        token = path.stem.rsplit("_", 1)[-1].upper()
        return token if token in rungs.RUNGS else None
    map_rung, states_rung = _named_rung(args.soil_map), _named_rung(args.states)
    if map_rung is not None and states_rung is not None and map_rung != states_rung:
        raise SystemExit(
            f"{args.soil_map.name} is {map_rung} and {args.states.name} is "
            f"{states_rung}. The states file is per build AND per rung and every "
            "consumer resolves it by rung, so writing one rung's rows under "
            "another's name gives ExoPlaSim and LPJ-GUESS a column count that is "
            "not their grid's. Name both at the same rung.")
    if not args.soil_map.is_file():
        raise SystemExit(
            f"{args.soil_map} does not exist. Run "
            "pedology/scripts/build_soil.py --grid "
            "source/<build>/exoplasim-<rung>, which weathers the lithology "
            "under a CLIMATOLOGY and so builds at whatever rung one is "
            "DECLARED at. This script derives from that soil map alone and "
            "adds no climate; see exoplasim/notes/route-step-criteria.md for "
            "why the soil at a rung waits for a climate at that rung rather "
            "than being built on one remapped onto it.")

    decl = load()
    written: list[str] = []
    if args.update_declaration:
        measured, _ = check_thermal_against_states(
            decl, read_soil_map(args.soil_map), gravity)
        written = update_declaration(measured)
        decl = load()
    report = build_report(decl, config, args.soil_map, gravity,
                          None if args.no_emit else args.states)
    report["declaration_written"] = written
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

        case = report["adopted_case"]
        adopted = case["adopted"]
        print(f"\n  adopted central case, world-of6n: {adopted['name']}")
        for who in ("exoplasim_dwmax", "lpj_guess_capacity"):
            move = adopted[who]
            print(f"      {who:20s} "
                  f"{move['median_mm_before']:7.1f} -> "
                  f"{move['median_mm_after']:7.1f} mm median, "
                  f"x{move['ratio']['median']:.3f} "
                  f"({move['ratio']['p10']:.3f} to {move['ratio']['p90']:.3f})")
        for name in case["not_adopted"]:
            print(f"      not adopted: {name}")

        frames = report["frame_consistency"]
        ratios = frames["available_ratio_to_shipped"]
        print(f"\n  frame consistency at {frames['gravity_ratio_to_earth']:.4f} "
              "Earth gravities, available capacity against the shipped code")
        for name in ("earth_equivalent_frame", "vesper_frame",
                     "frame_mixed_wilting_only"):
            print(f"    {name:26s} {ratios[name]['median']:.4f} median, "
                  f"{ratios[name]['p10']:.4f} to {ratios[name]['p90']:.4f}")
        print(f"    the two consistent frames agree to "
              f"{frames['frames_agree_max_abs_volumetric']:.1e} volumetric, and "
              "the frame-mixed variant does not agree with them")

        emitted = report["emitted_states"]
        if emitted is not None:
            print(f"\n  emitted {emitted['path']}")
            print(f"    awc                   {emitted['awc_mm']['median']:7.1f} mm "
                  f"median, {emitted['awc_mm']['p10']:.1f} to "
                  f"{emitted['awc_mm']['p90']:.1f} p10-p90 over "
                  f"{emitted['cells']} cells")

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

        unc = report["uncertainty_cases"]
        awc = unc["awc_mm"]
        print(f"\n  uncertainty of the adopted case, from cosby1984 Table 4's "
              f"within-class S.D. rows")
        print(f"    awc envelope          {awc['low']['median']:7.1f} to "
              f"{awc['high']['median']:7.1f} mm median, "
              f"x{awc['low_over_central']['median']:.3f} to "
              f"x{awc['high_over_central']['median']:.3f} of the central case")
        cmp_g = unc["against_the_gravity_shift"]
        print(f"    against the gravity shift  "
              f"{cmp_g['case_move_in_awc']['median']:.3f} against "
              f"{cmp_g['gravity_move_in_awc']['median']:.3f}, a factor "
              f"{cmp_g['case_over_gravity']['median']:.1f} at the median cell; "
              f"smaller on {cmp_g['cells_where_the_case_move_is_smaller']} cells")
        print(f"    ordering violations   "
              f"{unc['ordering_violations_over_all_cases']} over "
              f"{unc['corners']} cases and every cell")

        th = report["thermal_reduction"]
        ti = th["thermal_inertia_j_m2_k_s05"]
        print(f"\n  the declared thermal reduction, on this build's own column")
        print(f"    saturation the store spans "
              f"{th['sr_at_wilting_point']:.4f} empty to "
              f"{th['sr_at_field_capacity']:.4f} full")
        print(f"    thermal inertia            {ti['empty_store']:7.1f} empty to "
              f"{ti['full_store']:7.1f} full, a factor "
              f"{ti['store_range_factor']:.3f}")
        print(f"    the full saturation sweep  {ti['air_dry']:7.1f} air dry to "
              f"{ti['saturated']:7.1f} saturated, a factor "
              f"{ti['full_sweep_factor']:.3f}")
        print(f"    an empty store is the WILTING POINT, not a dry soil, so the "
              f"store cannot reach the dry end")

        if written:
            print(f"\n  wrote the derived saturation mapping into "
                  f"{DECLARATION.relative_to(PROJECT_ROOT)}")
            for line in written:
                print(f"    {line}")
        print(f"\n  the declared mapping is restated in {len(report['restatements'])} "
              "places, each with its own enforcer")
        for row in report["restatements"]:
            states = "-" if row["states"] is None else f"{row['states']:.4f}"
            print(f"    {row['status']:<7} {row['where']}  {row['key']} "
                  f"= {states} against {row['declared']:.4f}")
        stale = [row for row in report["restatements"]
                 if row["status"] == "stale"]
        if stale:
            print("    a stale restatement is refused by the enforcer named "
                  "against it, not here; update it in the same commit")
        if any(row["status"] == "unset" for row in report["restatements"]):
            print("    unset is one copy instead of two: the namelist key is "
                  "written at landmod.f90's compiled default, which is a row "
                  "above and is checked against this contract")

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
