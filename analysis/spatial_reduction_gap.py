#!/usr/bin/env python3
"""How much does the ORDER of a mesh-to-grid reduction change the answer? SPAT-7.

    python analysis/spatial_reduction_gap.py                 # every rung
    python analysis/spatial_reduction_gap.py --rung T42      # one of them

Worldbuilding. Vesper is an invented super-Earth; every quantity here is a
modelled field of that planet, measured on a World Orogen export. Nothing is
run: no climate model, no soil solve, no biosphere. This is arithmetic on
terrain and on the parameterisations the pipeline already applies to it.

## The question

`biosphere/notes/spatial-support-ecological-aggregation-audit.md` finding 4
lists paths where a cell mean is handed to a nonlinear law. Rastetter et al.
(1992) is the general result: `f(mean(x))` is not `mean(f(x))`, and the gap
between them is set by the curvature of `f` and the within-cell spread of `x`.
The audit names the paths and does not size them, so it cannot say which ones
are worth an operator. This sizes them.

The procedure is `hydrography/notes/subgrid-water-table.md` section 2, which a
sibling batch wrote for exactly this class and which is not re-invented here:
expand the parameterisation about the cell mean, COMPUTE the correction rather
than argue it, and compare it against the size of the effect the term carries.
Its second clause binds every arm below: sub-grid information may reach a
cell-scale parameter only as a statistic of the cell's own distribution -- a
fraction, a rank, a moment -- never as a resolved gradient or a within-cell
position.

## What is measured, and in whose units

Each arm reports the gap in the units the CONSUMING step reports, so a reader
can tell which reductions change a decision. Two arms are not Jensen gaps at
all and say so: they measure a calibration or a statistic that moves with the
support rather than with the operator order.

## The bars, declared before anything was run

A gap is MATERIAL when it exceeds the instrument the consuming step already
reports itself against, and NOISE when it does not. Each arm names its own
instrument, because a single percentage across six different quantities would
be a preference rather than a criterion:

  roughness_exchange      the step's own `z_ref` bracket. `build_surface_
                          roughness.py` reports every exchange coefficient as
                          a range over 273.15-313.15 K because the reference
                          height is not known before a climatology exists. A
                          gap inside that range cannot be distinguished from
                          the ignorance the step already declares.
  clay_conversion         a factor 1.35 in the weathering intensity, which is
                          Dunne's own S_y.x = 0.13 log units carried in
                          `pedology/config/pedogenesis.yaml` as the scatter on
                          the runoff exponent. A texture shift smaller than
                          what that scatter produces is inside the law's
                          published error.
  regolith_erodibility    0.02 m, `regolith.minimum_depth_m`: the thinnest
                          profile the pedogenesis model distinguishes at all.
  saturation_deficit      1% of the saturation vapour pressure. Declared, not
                          sourced: no estimator disagreement for this world's
                          evaporation is available at this step, and 1% is
                          below every water-balance term the carve verdict
                          weighs. Bracketed over lapse rate and reference
                          temperature and reported as a bracket.
  roughness_recalibration NOT a Jensen gap. The solved orographic coefficient
                          is a calibration that changes with the support, and
                          what is reported is how much of a cross-rung
                          difference belongs to it rather than to the terrain.
  subgrid_slope_support   NOT a Jensen gap. `pedology/scripts/build_soil.py:
                          subgrid_slope` divides a within-cell elevation
                          spread by the MESH spacing, so it carries a length
                          from the mesh; the arm measures the same cell on this
                          project's two builds and reports the ratio.

## What this cannot measure here

Every arm that needs a climate field is bracketed rather than evaluated: this
project's `baseline_climatology` is null and `hydrography/data/<build>/
water_table.nc` does not exist, so runoff, temperature and water table depth
have no values to read. Where an arm needs one it sweeps a declared range and
reports the whole range. An arm whose gap is material across the WHOLE bracket
is material whatever the climate turns out to be; one that crosses its bar
inside the bracket is reported as crossing and is not resolved here.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "lib"))
sys.path.insert(0, str(PROJECT_ROOT / "exoplasim" / "scripts"))

import builds                                                  # noqa: E402
from gridding import RUNGS, region_cells                       # noqa: E402
from orogen import Export, LAND                                # noqa: E402
from lapse import reference_height_m                           # noqa: E402
# ONE definition of how much forest the land carries, and one of the height a
# bulk transfer coefficient is taken over. Both are imported from the step that
# owns them rather than restated, which is `failure-modes.md` class 17: a second
# formulation of one quantity is not a check on the first.
from build_surface_roughness import (CE_BRACKET_K, DEFAULT_BARE_Z0_M,   # noqa: E402
                                     DEFAULT_CANOPY_Z0_M,
                                     DEFAULT_FOREST_Z0_M,
                                     EXOPLASIM_DZ0LAND_M, KARMAN)
from build_surface_albedo import MODE_FOREST_FRACTION           # noqa: E402

OUTPUT = PROJECT_ROOT / "analysis" / "spatial_reduction_gap.json"

# The intensity bracket the weathering arm sweeps: the clip range
# `pedogenesis.yaml` declares for W itself, so the sweep covers every value the
# law can return rather than a guess at this world's climate.
INTENSITY_SWEEP = (0.02, 0.2, 1.0, 2.0, 6.0)
DUNNE_SCATTER = 1.35        # S_y.x = 0.13 log units, pedogenesis.yaml

# The regolith arm's one free ratio, `erosion_weight * E / P`, swept over the
# decades the saturating depth law can reach: at 0.01 the profile sits at its
# ceiling and at 100 it is scraped to bare rock, so the sweep spans the whole
# range of the law rather than a guess at this world's climate.
RHO_SWEEP = (0.01, 0.1, 0.5, 1.0, 3.0, 10.0, 100.0)

# Lapse rate bracket for the saturation-deficit arm, K per km. `lib/lapse.py`
# owns the rates this world runs at; the bracket spans dry-adiabatic-ish to
# strongly moist, because which one applies is a property of a climatology this
# step does not have.
LAPSE_BRACKET_K_PER_KM = (4.0, 9.8)
REFERENCE_AIR_BRACKET_K = CE_BRACKET_K       # the same liquid-water span
SATURATION_BAR = 0.01                        # 1% of e_sat, declared above

L_VAP = 2.5e6               # J/kg, latent heat of vaporisation
R_VAP = 461.5               # J/kg/K, gas constant for water vapour


# --------------------------------------------------------------------------
# arm 1: roughness, and the exchange coefficient it acts through
# --------------------------------------------------------------------------

def _ce(z0, z_ref: float) -> np.ndarray:
    """Neutral bulk exchange coefficient. `build_surface_roughness.py`'s own."""
    return KARMAN ** 2 / np.log(z_ref / np.maximum(np.asarray(z0), 1e-6)) ** 2


def _solve_k_oro(z0_surf_cell, sigma_m, weight, target: float) -> float:
    """The orographic coefficient `build_surface_roughness.py` solves for.

    Reproduced here rather than imported because the builder solves it inside
    `main()`; the bisection is the builder's, on the builder's own quantities.
    """
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        trial = float((np.sqrt(z0_surf_cell ** 2 + (mid * sigma_m) ** 2)
                       * weight).sum() / weight.sum())
        if trial < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def roughness_arm(mesh: Export, cell, ncell, land, area, config) -> dict:
    """Aggregate-then-process against process-then-aggregate on `ce`.

    The builder averages `z0_surface` over a cell's land, adds the cell's
    orographic term in quadrature, and hands ONE roughness to the model, which
    then takes `ce = k^2 / ln(z_ref/z0)^2`. `ce` is what the turbulent flux is
    linear in, so a relative gap in `ce` is a relative gap in the evaporation
    and the sensible heat the cell delivers.

    The process-then-aggregate arm gives every mesh region its own surface
    roughness, adds the SAME cell orographic term -- subgrid relief is a cell
    statistic and has no per-region value, which is the sub-grid rule's second
    clause -- takes `ce` per region, and area-averages that.
    """
    model = config["model"]
    rock = mesh.field("substrate_class").astype(int)
    codes = [r["code"] for r in mesh.manifest["lithology"]["rockClasses"]]
    barren = np.zeros(rock.shape, dtype=bool)
    for code in model.get("barren_rock_classes", []):
        if code in codes:
            barren |= rock == codes.index(code)

    fraction = MODE_FOREST_FRACTION.get(str(model.get("land_albedo_source",
                                                      "lithology")))
    canopy = DEFAULT_CANOPY_Z0_M
    if fraction:
        canopy = (1.0 - fraction) * DEFAULT_CANOPY_Z0_M + fraction * DEFAULT_FOREST_Z0_M
    z0_surface = np.where(barren, DEFAULT_BARE_Z0_M, canopy)

    elev_km = mesh.field("elevation_km").astype(np.float64)
    w = np.bincount(cell[land], weights=area[land], minlength=ncell)
    bare = np.bincount(cell[land & barren], weights=area[land & barren],
                       minlength=ncell)
    m1 = np.bincount(cell[land], weights=area[land] * z0_surface[land], minlength=ncell)
    e1 = np.bincount(cell[land], weights=area[land] * elev_km[land], minlength=ncell)
    e2 = np.bincount(cell[land], weights=area[land] * elev_km[land] ** 2, minlength=ncell)
    have = w > 0
    z0_surf_cell = np.where(have, m1 / np.maximum(w, 1e-30), 0.0)
    emean = np.where(have, e1 / np.maximum(w, 1e-30), 0.0)
    sigma_m = np.sqrt(np.maximum(np.where(have, e2 / np.maximum(w, 1e-30)
                                          - emean ** 2, 0.0), 0.0)) * 1000.0

    k_oro = _solve_k_oro(z0_surf_cell[have], sigma_m[have], w[have],
                         EXOPLASIM_DZ0LAND_M)

    z0_cell = np.sqrt(z0_surf_cell ** 2 + (k_oro * sigma_m) ** 2)
    z0_region = np.sqrt(z0_surface[land] ** 2 + (k_oro * sigma_m[cell[land]]) ** 2)

    # The instrument, computed FIRST because every arm below is judged against
    # it: the step's own bracket on `ce` at the land-mean roughness, which is
    # the range `build_surface_roughness.py` already reports because `z_ref` is
    # not known before a climatology exists.
    ce_bracket = [float(_ce(np.array([np.average(z0_cell[have], weights=w[have])]),
                            reference_height_m(t, config))[0])
                  for t in CE_BRACKET_K]
    instrument = abs(ce_bracket[1] - ce_bracket[0]) / min(ce_bracket)

    out = {"orographic_coefficient_solved": round(k_oro, 6),
           "land_cells": int(have.sum()),
           "subgrid_stdev_m_land_median": round(float(np.median(sigma_m[have])), 2),
           "barren_land_area_fraction": round(
               float(area[land & barren].sum() / area[land].sum()), 6),
           "z_ref_bracket_m": [round(reference_height_m(t, config), 2)
                               for t in CE_BRACKET_K],
           "by_reference_air_k": {}}

    for t_air in CE_BRACKET_K:
        z_ref = reference_height_m(t_air, config)
        ce_of_mean = _ce(z0_cell, z_ref)
        ce_region = _ce(z0_region, z_ref)
        num = np.bincount(cell[land], weights=area[land] * ce_region, minlength=ncell)
        mean_of_ce = np.where(have, num / np.maximum(w, 1e-30), 0.0)

        gap = mean_of_ce - ce_of_mean
        land_area = w[have]
        agg = float((ce_of_mean[have] * land_area).sum() / land_area.sum())
        exp = float((mean_of_ce[have] * land_area).sum() / land_area.sum())
        rel = np.where(ce_of_mean[have] > 0, gap[have] / ce_of_mean[have], 0.0)
        # STRATIFIED BY BARREN SHARE, and the bar does not move. The tail was
        # found first and this is the explanation being tested, not a second
        # criterion: the barren classes carry a roughness two orders below the
        # canopy, so a cell that is mostly barren has its exchange coefficient
        # set by the minority surface once the roughnesses are averaged. Those
        # are the closed-basin floors the carve verdict integrates evaporation
        # over, which is why the roughness field exists at all.
        bare_share = np.where(have, bare / np.maximum(w, 1e-30), 0.0)[have]
        strata = {}
        for lo_q, hi_q in ((0.0, 0.5), (0.5, 0.9), (0.9, 1.0)):
            lo_v, hi_v = np.quantile(bare_share, [lo_q, hi_q])
            sel = (bare_share >= lo_v) & (bare_share <= hi_v)
            if not sel.any():
                continue
            wa = land_area[sel]
            strata[f"barren_share_q{lo_q}_{hi_q}"] = {
                "barren_share_range": [round(float(lo_v), 4), round(float(hi_v), 4)],
                "land_area_weighted_relative_gap": round(
                    float((rel[sel] * wa).sum() / wa.sum()), 6),
                "per_cell_relative_gap": {
                    str(p): round(float(np.percentile(rel[sel], p)), 6)
                    for p in (1, 50, 99)},
            }
        out["by_reference_air_k"][str(t_air)] = {
            "ce_aggregate_then_process_land_mean": round(agg, 6),
            "ce_process_then_aggregate_land_mean": round(exp, 6),
            "land_mean_relative_gap": round((exp - agg) / agg, 6),
            "per_cell_relative_gap": {
                str(p): round(float(np.percentile(rel, p)), 6)
                for p in (1, 25, 50, 75, 95, 99)},
            "land_area_fraction_gap_beyond_instrument": round(
                float(land_area[np.abs(rel) > instrument].sum() / land_area.sum()), 6),
            "by_barren_share": strata,
        }

    out["instrument_z_ref_bracket_relative"] = round(instrument, 6)
    return out


def recalibration_arm(per_rung: dict, reference: str) -> dict:
    """What the per-rung solve of the orographic coefficient is worth.

    `build_surface_roughness.py` re-solves the coefficient at every resolution
    so that the land mean lands on `dz0land`. That is a calibration whose value
    depends on the support, so a difference between two rungs mixes the terrain
    with the calibration and finding 8 of the audit cannot be answered from it.
    This reports the coefficient per rung and the land-mean roughness each rung
    would carry if the reference rung's coefficient were held fixed instead.
    """
    ref = per_rung[reference]["orographic_coefficient_solved"]
    rows = {}
    for rung, arm in per_rung.items():
        k = arm["orographic_coefficient_solved"]
        rows[rung] = {
            "coefficient_solved": k,
            "ratio_to_reference": round(k / ref, 6) if ref else None,
            "subgrid_stdev_m_land_median": arm["subgrid_stdev_m_land_median"],
        }
    return {"reference_rung": reference, "target_land_mean_m": EXOPLASIM_DZ0LAND_M,
            "by_rung": rows}


# --------------------------------------------------------------------------
# arm 2: mixing parent rock before weathering it
# --------------------------------------------------------------------------

def clay_arm(mesh: Export, cell, ncell, land, area, texture_cfg) -> dict:
    """Mix-then-weather against weather-then-mix, in texture mass fractions.

    `pedology/scripts/build_soil.py:weather_texture` area-mixes every rock
    class's parent texture inside a cell and THEN converts weatherable minerals
    to clay. All three fractions are measured, not clay alone, because the
    operator is not equally nonlinear in them: the clay TOTAL is
    `clay + clip(1 - quartz - clay) * yield * sat(W)`, affine in the mixed
    quartz and clay shares wherever the clip does not bind, while the SPLIT of
    that conversion between sand and silt goes through a ratio and two `minimum`
    clamps and is not. So this arm reproduces `weather_texture`'s arithmetic
    exactly rather than its clay line only.

    The intensity is swept over its own declared clip range rather than read
    from a climatology, because this world has none yet. Texture is what the
    soil map reports and what plant-available capacity is built from.
    """
    rock = mesh.field("substrate_class").astype(int)
    entries = mesh.manifest["lithology"]["rockClasses"]
    yield_ = float(texture_cfg["clay_yield"])
    conv = float(texture_cfg["clay_conversion"])

    # Per-region parent texture, from the region's own rock class.
    n_class = len(entries)
    parent = np.zeros((n_class, 4))
    for i, entry in enumerate(entries):
        t = texture_cfg.get(entry["code"])
        if t is None:
            raise SystemExit(
                f"rock class {entry['code']!r} has no texture in "
                "pedogenesis.yaml; this arm cannot be measured without it")
        parent[i] = (t["sand"], t["silt"], t["clay"], t["quartz"])

    idx = rock[land]
    a = area[land]
    c = cell[land]
    w = np.bincount(c, weights=a, minlength=ncell)
    have = w > 0

    # Cell-mean parent texture: the builder's own mix.
    mixed = np.zeros((ncell, 4))
    for k in range(4):
        mixed[:, k] = np.bincount(c, weights=a * parent[idx, k], minlength=ncell)
    mixed[have] /= w[have, None]

    ratio = float(texture_cfg["sand_to_silt_loss_ratio"])

    def weathered(sand, silt, clay, quartz, saturating):
        """`weather_texture`'s arithmetic, line for line, on any support."""
        weatherable = np.clip(1.0 - quartz - clay, 0.0, 1.0) * yield_
        converted = weatherable * saturating
        denominator = np.maximum(sand * ratio + silt, 1e-9)
        from_sand = np.minimum(converted * sand * ratio / denominator, sand)
        from_silt = np.minimum(converted * silt / denominator, silt)
        return (np.clip(sand - from_sand, 0.0, 1.0),
                np.clip(silt - from_silt, 0.0, 1.0),
                np.clip(clay + from_sand + from_silt, 0.0, 1.0))

    rows = {}
    land_area = w[have]
    for intensity in INTENSITY_SWEEP:
        saturating = 1.0 - np.exp(-conv * intensity)
        agg = weathered(mixed[:, 0], mixed[:, 1], mixed[:, 2], mixed[:, 3],
                        saturating)
        per_region = weathered(parent[idx, 0], parent[idx, 1], parent[idx, 2],
                               parent[idx, 3], saturating)
        row = {}
        for k, name in enumerate(("sand", "silt", "clay")):
            exp = np.bincount(c, weights=a * per_region[k], minlength=ncell)
            exp[have] /= w[have]
            gap = exp[have] - agg[k][have]
            row[name] = {
                "aggregate_then_process_land_mean": round(
                    float((agg[k][have] * land_area).sum() / land_area.sum()), 6),
                "process_then_aggregate_land_mean": round(
                    float((exp[have] * land_area).sum() / land_area.sum()), 6),
                "land_mean_gap": round(
                    float((gap * land_area).sum() / land_area.sum()), 8),
                "per_cell_gap": {str(p): round(float(np.percentile(gap, p)), 8)
                                 for p in (1, 50, 99)},
                "max_abs_cell_gap": round(float(np.abs(gap).max()), 8),
            }
        rows[str(intensity)] = row

    # The instrument: what Dunne's own scatter on the runoff exponent does to
    # the same cell's texture, at each swept intensity.
    instrument = {}
    for intensity in INTENSITY_SWEEP:
        base = weathered(mixed[have, 0], mixed[have, 1], mixed[have, 2],
                         mixed[have, 3], 1.0 - np.exp(-conv * intensity))
        hi = weathered(mixed[have, 0], mixed[have, 1], mixed[have, 2],
                       mixed[have, 3],
                       1.0 - np.exp(-conv * intensity * DUNNE_SCATTER))
        instrument[str(intensity)] = {
            name: round(float((np.abs(hi[k] - base[k]) * land_area).sum()
                              / land_area.sum()), 8)
            for k, name in enumerate(("sand", "silt", "clay"))}
    return {"intensity_sweep": list(INTENSITY_SWEEP),
            "clay_yield": yield_, "clay_conversion": conv,
            "sand_to_silt_loss_ratio": ratio,
            "by_intensity": rows,
            "instrument_dunne_scatter": DUNNE_SCATTER,
            "instrument_land_mean_shift": instrument}


# --------------------------------------------------------------------------
# arm 3: mixing erodibility before the saturating depth law
# --------------------------------------------------------------------------

def regolith_arm(mesh: Export, cell, ncell, land, area, pedo) -> dict:
    """Mixing erodibility before the saturating depth law, in metres of regolith.

    `build_soil.py` area-mixes every rock class's erodibility inside a cell and
    then applies

        depth = maximum_depth * P / (P + erosion_weight * E),  E ~ erodibility

    which is convex in the erodibility, so the mean of the depths is above the
    depth at the mean. Erodibility is a per-region property of the rock the
    region is made of, so both orders are well defined and the comparison is
    the operator's and nothing else's.

    Everything the law needs except the erodibility collapses into ONE ratio,

        rho = erosion_weight * E(cell mean erodibility) / P

    so no climatology is needed to size this: the arm sweeps `rho` over the
    decades the law can reach and reports the gap in metres at each. `rho` is
    the fraction of the way the cell has been driven from the ceiling towards
    bare rock; `depth = maximum_depth / (1 + rho)` at the cell mean.
    """
    d_min = float(pedo["regolith"]["minimum_depth_m"])
    d_max = float(pedo["regolith"]["maximum_depth_m"])

    entries = mesh.manifest["lithology"]["rockClasses"]
    erod_of_class = np.array([float(e["erodibility"]) for e in entries])
    rock = mesh.field("substrate_class").astype(int)
    erod = erod_of_class[rock]

    a = area[land]
    c = cell[land]
    w = np.bincount(c, weights=a, minlength=ncell)
    have = w > 0
    mixed = np.bincount(c, weights=a * erod[land], minlength=ncell)
    mixed[have] /= w[have]

    rows = {}
    for rho in RHO_SWEEP:
        # The cell-mean arm, which is what the builder computes.
        agg = d_max / (1.0 + rho)
        # The per-region arm: every region carries its own rock's erodibility,
        # scaled so the cell mean is the same `rho`.
        with np.errstate(invalid="ignore", divide="ignore"):
            scaled = np.where(mixed[c] > 0, erod[land] / mixed[c], 1.0)
        per_region = d_max / (1.0 + rho * scaled)
        exp = np.bincount(c, weights=a * per_region, minlength=ncell)
        exp[have] /= w[have]
        gap = np.clip(exp[have], d_min, d_max) - np.clip(agg, d_min, d_max)
        land_area = w[have]
        rows[str(rho)] = {
            "depth_aggregate_then_process_m": round(float(np.clip(agg, d_min, d_max)), 6),
            "depth_process_then_aggregate_land_mean_m": round(
                float((np.clip(exp[have], d_min, d_max) * land_area).sum()
                      / land_area.sum()), 6),
            "land_mean_gap_m": round(float((gap * land_area).sum() / land_area.sum()), 6),
            "per_cell_gap_m": {str(p): round(float(np.percentile(gap, p)), 6)
                               for p in (50, 75, 95, 99)},
            "land_area_fraction_above_bar": round(
                float(land_area[np.abs(gap) > d_min].sum() / land_area.sum()), 6),
        }
    return {"bar_m": d_min, "maximum_depth_m": d_max,
            "erodibility_class_range": [round(float(erod_of_class.min()), 3),
                                        round(float(erod_of_class.max()), 3)],
            "land_mean_erodibility": round(
                float((erod[land] * a).sum() / a.sum()), 4),
            "rho_sweep": list(RHO_SWEEP), "by_rho": rows}


# --------------------------------------------------------------------------
# arm 4: subgrid elevation through Clausius-Clapeyron
# --------------------------------------------------------------------------

def saturation_arm(mesh: Export, cell, ncell, land, area) -> dict:
    """What a mean-elevation cell costs in saturation vapour pressure.

    `build_boundary_conditions.py` averages elevation over a cell's land and
    hands the model one height. Temperature is LINEAR in height through a lapse
    rate, so that reduction is exact for temperature and the Jensen term is
    zero. Saturation vapour pressure is not: it is exponential in temperature,
    so the cell's own spread of heights raises the mean saturation deficit
    above the one computed at the mean height. That is the term a single
    elevation per cell discards, and evaporation is what consumes it.

    Bracketed twice over, and both brackets are reported: the lapse rate,
    because which one applies is a property of a climatology this step does not
    have, and the reference air temperature, over the span in which surface
    water is liquid.
    """
    elev_km = mesh.field("elevation_km").astype(np.float64)
    a = area[land]
    c = cell[land]
    w = np.bincount(c, weights=a, minlength=ncell)
    have = w > 0
    m1 = np.bincount(c, weights=a * elev_km[land], minlength=ncell)
    mean = np.where(have, m1 / np.maximum(w, 1e-30), 0.0)

    rows = {}
    for lapse in LAPSE_BRACKET_K_PER_KM:
        for t0 in REFERENCE_AIR_BRACKET_K:
            # e_sat ratio against the cell mean height, per region.
            dt = -lapse * (elev_km[land] - mean[c])
            ratio = np.exp(L_VAP / R_VAP * (1.0 / t0 - 1.0 / (t0 + dt)))
            num = np.bincount(c, weights=a * ratio, minlength=ncell)
            gap = np.where(have, num / np.maximum(w, 1e-30), 1.0) - 1.0
            land_area = w[have]
            rows[f"lapse{lapse}_t{t0}"] = {
                "land_mean_relative_gap": round(
                    float((gap[have] * land_area).sum() / land_area.sum()), 6),
                "per_cell_relative_gap": {
                    str(p): round(float(np.percentile(gap[have], p)), 6)
                    for p in (50, 75, 95, 99)},
                "land_area_fraction_above_bar": round(
                    float(land_area[gap[have] > SATURATION_BAR].sum()
                          / land_area.sum()), 6),
            }
    return {"bar_relative": SATURATION_BAR,
            "lapse_bracket_k_per_km": list(LAPSE_BRACKET_K_PER_KM),
            "reference_air_bracket_k": list(REFERENCE_AIR_BRACKET_K),
            "by_bracket_corner": rows}


# --------------------------------------------------------------------------
# arm 5: the support a length-carrying statistic is measured on
# --------------------------------------------------------------------------

def support_arm(config, rung: str, pedo) -> dict:
    """`subgrid_slope` on this project's two builds, at the same grid.

    The statistic is an elevation spread over the MESH spacing, so it carries a
    length from the mesh rather than from the terrain. Both builds are the same
    planet at the same seed; the ratio between them is what changing the region
    count alone does to a quantity the pedogenesis model consumes as a
    gradient. Reported as the ratio and as what it is worth in metres of
    regolith through the catena divisor, which is the consuming step.
    """
    s = float(pedo["catena"]["slope_transport"])
    out = {}
    for name in ("precarve-craton", "precarve-craton-10m"):
        root = PROJECT_ROOT / "source" / name
        mesh_dir = builds.mesh_export_of(root)
        grid_dir = root / f"exoplasim-{rung}"
        if not mesh_dir.is_dir() or not grid_dir.is_dir():
            out[name] = {"absent": True}
            continue
        mesh = Export(mesh_dir)
        cell, nlat, nlon = region_cells(mesh, grid_dir)
        ncell = nlat * nlon
        land = mesh.surface_class == LAND
        elev_m = mesh.field("elevation_km").astype(np.float64) * 1000.0
        spacing_m = float(mesh.manifest["basins"]["resolution"]["avgEdgeKm"]) * 1000.0
        c = cell[land]
        n = np.bincount(c, minlength=ncell).astype(np.float64)
        t1 = np.bincount(c, weights=elev_m[land], minlength=ncell)
        t2 = np.bincount(c, weights=elev_m[land] ** 2, minlength=ncell)
        with np.errstate(invalid="ignore", divide="ignore"):
            m = np.where(n > 0, t1 / np.maximum(n, 1), 0.0)
            v = np.where(n > 1, t2 / np.maximum(n, 1) - m ** 2, 0.0)
        tanb = np.sqrt(np.maximum(v, 0.0)) / spacing_m
        have = n > 1
        out[name] = {
            "regions": int(mesh.n_regions),
            "mesh_spacing_km": round(spacing_m / 1000.0, 4),
            "land_median_tan_beta": round(float(np.median(tanb[have])), 6),
            "land_median_stdev_m": round(float(np.median(np.sqrt(
                np.maximum(v[have], 0.0)))), 3),
            "land_median_catena_divisor": round(
                float(np.median(1.0 / (1.0 + s * tanb[have]))), 6),
        }
    if all("absent" not in v for v in out.values()):
        a, b = out["precarve-craton"], out["precarve-craton-10m"]
        out["ratio_10m_over_2p5m"] = {
            "tan_beta": round(b["land_median_tan_beta"] / a["land_median_tan_beta"], 4),
            "stdev_m": round(b["land_median_stdev_m"] / a["land_median_stdev_m"], 4),
            "catena_divisor": round(
                b["land_median_catena_divisor"] / a["land_median_catena_divisor"], 4),
        }
    return {"rung": rung, "slope_transport": s, "builds": out}


# --------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path,
                    default=PROJECT_ROOT / "config" / "planet.yaml")
    ap.add_argument("--rung", action="append", default=None,
                    help="ladder rung to measure; repeatable, defaults to "
                         "every rung the configured build carries an export for")
    ap.add_argument("--reference-rung", default="T42",
                    help="rung the roughness recalibration is compared against")
    ap.add_argument("--output", type=Path, default=OUTPUT)
    args = ap.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    pedo = yaml.safe_load(
        (PROJECT_ROOT / "pedology" / "config" / "pedogenesis.yaml")
        .read_text(encoding="utf-8"))

    build_root = builds.build_root(config)
    mesh = Export(builds.mesh_export(config))
    area = mesh.cell_area.astype(np.float64)
    land = mesh.surface_class == LAND

    rungs = args.rung or [r for r in RUNGS
                          if (build_root / f"exoplasim-{r}").is_dir()]
    per_rung: dict[str, dict] = {}
    for rung in rungs:
        grid_dir = build_root / f"exoplasim-{rung}"
        cell, nlat, nlon = region_cells(mesh, grid_dir)
        ncell = nlat * nlon
        per_rung[rung] = {
            "grid": [nlat, nlon],
            "roughness_exchange": roughness_arm(mesh, cell, ncell, land, area, config),
            "clay_conversion": clay_arm(mesh, cell, ncell, land, area,
                                        pedo["texture"]),
            "regolith_erodibility": regolith_arm(mesh, cell, ncell, land,
                                                area, pedo),
            "saturation_deficit": saturation_arm(mesh, cell, ncell, land, area),
        }
        print(f"{rung}: done")

    reference = (args.reference_rung if args.reference_rung in per_rung
                 else sorted(per_rung)[0])
    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "issue": "SPAT-7",
        "build": build_root.name,
        "terrain_hash": mesh.terrain_hash,
        "regions": int(mesh.n_regions),
        "bars": {
            "roughness_exchange": "the step's own z_ref bracket over "
                                  "273.15-313.15 K",
            "clay_conversion": f"a factor {DUNNE_SCATTER} in the weathering "
                               "intensity, Dunne's S_y.x carried in "
                               "pedogenesis.yaml",
            "regolith_erodibility": "regolith.minimum_depth_m",
            "saturation_deficit": f"{SATURATION_BAR} of e_sat, declared",
        },
        "by_rung": {k: {n: v for n, v in arm.items()}
                    for k, arm in per_rung.items()},
        "roughness_recalibration": recalibration_arm(
            {k: v["roughness_exchange"] for k, v in per_rung.items()}, reference),
        "subgrid_slope_support": support_arm(config, reference, pedo),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
