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
                          temperature and reported as a bracket. The bar has
                          NOT been moved since the first measurement crossed
                          it: there is no sourced replacement to move it to,
                          and a bar placed after the result it judges is not a
                          criterion. What was corrected is the bracket, and
                          what the corrections are is fixed below.
  roughness_recalibration NOT a Jensen gap. The orographic term is DERIVED and
                          carries no solved coefficient, so a cross-rung
                          difference in the land-mean roughness is terrain and
                          aggregation alone; the arm reports the spread across
                          the ladder, which is the claim that can fail.
  subgrid_slope_support   NOT a Jensen gap. A within-cell elevation spread
                          quoted over the MESH spacing carries a length from
                          the mesh; the arm measures the same cell on this
                          project's two builds and reports the ratio over that
                          run and over the declared one.

## The saturation arm's corrections, and the verdict rule, fixed in advance

Three things in the arm were wrong for this planet or for this model, and each
is a defect independent of which way it moves the answer. Two of them are
corrected here; the third is reported rather than repaired because nothing at
this step can repair it.

1. THE LAPSE CEILING WAS EARTH'S. The bracket ran to 9.8 K per km, which is
   Earth's `g/cp` and not this planet's. `lib/lapse.py` owns the quantity and
   already treats the configured dry adiabat as the ceiling a measured rate
   must fall below, raising when it does not, so that is the ceiling the
   bracket takes. This planet's is steeper, so the correction WIDENS the
   bracket and moves the arm further from a verdict rather than closer to one.
2. THE FUNCTION WAS NOT THE MODEL'S. The arm evaluated an idealised
   Clausius-Clapeyron with a constant latent heat. The model evaluates
   Magnus-Teten with two coefficient sets and a phase switch at `tmelt`
   (`plasimmod.f90:ra1s/ra2s/ra4s`), and at the cold corner of the reference
   bracket a cell's high ground falls on the ice side, where the coefficients
   are steeper. The arm now evaluates the model's own function, with the
   coefficients READ from `p_earth.f90` in the manner `lib/sea_water.py`
   established rather than copied.
3. THE FLOOR HAS NOTHING UNDER IT. 4.0 K per km is a declared floor, not a
   bound: `lib/lapse.py` states that the measured environmental rate sits
   BELOW the moist adiabatic rate evaluated at the window-mean state, which is
   why it deliberately declines to floor there. Nothing at this step bounds
   the land-mean lapse rate away from zero, and the gap goes to zero with it.
   The low end of this arm's bracket is therefore a choice and the arm reports
   it as one.

The verdict rule, fixed before the corrected arm was run:

  - at or below the bar at EVERY corner of the corrected bracket, at some
    rung: admissible at that rung and above, no operator needed;
  - above the bar at EVERY corner: material whatever the climate turns out to
    be, and the correction is a sub-grid orographic term in the surface
    evaporation rather than a finer rung;
  - crossing the bar inside the bracket: NO VERDICT. What settles it is
    `lapse.environmental_lapse_k_per_km` on an accepted baseline climatology,
    which is the one input the corner sweep is standing in for. It is not
    settled by a finer rung, because refinement moves the gap by less than the
    bracket does.

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
import re
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
from lapse import dry_adiabat_k_per_km, reference_height_m     # noqa: E402
# ONE definition of how much forest the land carries, and one of the height a
# bulk transfer coefficient is taken over. Both are imported from the step that
# owns them rather than restated, which is `failure-modes.md` class 17: a second
# formulation of one quantity is not a check on the first.
from build_surface_roughness import (CE_BRACKET_K, DEFAULT_BARE_Z0_M,   # noqa: E402
                                     DEFAULT_CANOPY_Z0_M,
                                     DEFAULT_FOREST_Z0_M,
                                     EXOPLASIM_DZ0LAND_M, KARMAN,
                                     RESOLVED_BAND_MULTIPLE,
                                     effective_length, inner_layer_depth,
                                     orographic_form_drag,
                                     pressure_scale_height)
from build_surface_albedo import (MODE_FOREST_FRACTION, ROCK_BANDS,  # noqa: E402
                                  WETTING, band_shapes, rock_band_ratios,
                                  rock_wetting_ratios, _rock_id)
# The reduction operators, the mixing law and the climatology door, each from
# the module that owns it. `sadeghi_mix` is imported rather than restated
# because the compiled model runs the same shape and `soil_albedo_wetting.py`
# already checks the two against each other.
from gridding import (cell_expectation, cell_mean, gaussian_latitudes,  # noqa: E402
                      require_same_rows)
from climatology import annual_mean                             # noqa: E402
from paths import best_available_climatology                    # noqa: E402
from stellar import band_fractions                              # noqa: E402
sys.path.insert(0, str(PROJECT_ROOT / "analysis"))
from soil_albedo_wetting import (kubelka_munk_transform,        # noqa: E402
                                 sadeghi_mix)

OUTPUT = PROJECT_ROOT / "analysis" / "spatial_reduction_gap.json"

# The wetting arm the staged pair is built with. `build_surface_albedo.py`
# offers two and the field is staged from one; measuring the reduction on a
# different arm from the one that is staged would compare two things at once.
WETTING_ARM = "twomey"

# A cell of one substrate class has no spread for the mixing to be nonlinear
# over, so its gap is exactly zero. The tolerance is float64 round-off on an
# albedo, which is what an exact identity leaves behind.
ONE_CLASS_TOLERANCE = 1.0e-12

# The accepted baseline's state-storage tolerance, which
# `config/partial_surface.yaml` selected the tile operator against and
# `ocean/config/transport_loop.yaml` carries as the transport loop's own
# controlling scalar. An albedo error reaches it as `delta_alpha * S_down`.
STORAGE_TOLERANCE_W_M2 = 0.12

# The intensity bracket the weathering arm sweeps: the clip range
# `pedogenesis.yaml` declares for W itself, so the sweep covers every value the
# law can return rather than a guess at this world's climate.
INTENSITY_SWEEP = (0.02, 0.2, 1.0, 2.0, 6.0)
# Dunne's own regression standard error about the runoff fit, S_y.x, in log
# units. Paper-sourced, so nothing on this tree computes it; `pedogenesis.yaml`
# quotes the same figure and its antilog. The ANTILOG is arithmetic and is done
# here rather than restated, so this file carries one statement of the source
# number and none of the conversion.
DUNNE_SCATTER_LOG_UNITS = 0.13     # S_y.x, pedogenesis.yaml
DUNNE_SCATTER = 10.0 ** DUNNE_SCATTER_LOG_UNITS

# The regolith arm's one free ratio, `E / P`, swept over the
# decades the saturating depth law can reach: at 0.01 the profile sits at its
# ceiling and at 100 it is scraped to bare rock, so the sweep spans the whole
# range of the law rather than a guess at this world's climate.
RHO_SWEEP = (0.01, 0.1, 0.5, 1.0, 3.0, 10.0, 100.0)

# Lapse rate bracket for the saturation-deficit arm, K per km. `lib/lapse.py`
# owns the rates this world runs at and the ceiling comes from it, because a
# stably stratified column cannot lapse faster than THIS planet's dry adiabat
# and `lapse.environmental_lapse_k_per_km` already raises on a measured rate
# that does. The floor is DECLARED and has nothing under it: lapse.py states
# that the measured rate sits below the moist adiabatic rate at the window-mean
# state, so no moist floor is available, and the gap goes to zero with the rate.
LAPSE_FLOOR_K_PER_KM = 4.0
REFERENCE_AIR_BRACKET_K = CE_BRACKET_K       # the same liquid-water span
SATURATION_BAR = 0.01                        # 1% of e_sat, declared above

# The model's saturation vapour pressure is Magnus-Teten with two coefficient
# sets and a phase switch at `tmelt`, not an idealised Clausius-Clapeyron. The
# coefficients are READ from the planet module rather than copied, which is what
# `lib/sea_water.py` established: a copied model constant goes stale silently
# when the model moves, and this set moved once already when the ice branch was
# added.
PLANET_SOURCE = (PROJECT_ROOT / "vendor" / "exoplasim" / "exoplasim" / "plasim"
                 / "src" / "p_earth.f90")
MAGNUS_TETEN_NAMES = ("ra1", "ra2", "ra4", "ra1i", "ra2i", "ra4i", "tmelt")


def magnus_teten() -> dict:
    """The model's Magnus-Teten coefficients, read from `p_earth.f90`."""
    text = PLANET_SOURCE.read_text(encoding="utf-8")
    out = {}
    for name in MAGNUS_TETEN_NAMES:
        m = re.search(rf"^\s*{name}\s*=\s*([0-9.eEdD+-]+)", text, re.MULTILINE)
        if m is None:
            raise RuntimeError(
                f"{PLANET_SOURCE.name} no longer declares {name}; the arm "
                "cannot evaluate the model's saturation function")
        out[name] = float(m.group(1).replace("d", "e").replace("D", "e"))
    return out


def model_e_sat(t_k: np.ndarray, c: dict) -> np.ndarray:
    """`plasimmod.f90:ra1s/ra2s/ra4s` and `ra4d`, in Pa, elementwise.

    The phase follows the temperature, as it does at every land saturation site
    (`landmod.f90:716,932`); the denominator carries the same floor `ra4d`
    applies. A sea surface does NOT select this way, and this arm is land only.
    """
    t = np.asarray(t_k, dtype=np.float64)
    ice = t < c["tmelt"]
    a1 = np.where(ice, c["ra1i"], c["ra1"])
    a2 = np.where(ice, c["ra2i"], c["ra2"])
    a4 = np.where(ice, c["ra4i"], c["ra4"])
    return a1 * np.exp(a2 * (t - c["tmelt"]) / np.maximum(t - a4, 1.0))


def lapse_bracket_k_per_km(config: dict) -> tuple[float, float]:
    """Declared floor to this planet's own dry adiabat."""
    return (LAPSE_FLOOR_K_PER_KM, dry_adiabat_k_per_km(config))


# --------------------------------------------------------------------------
# arm 1: roughness, and the exchange coefficient it acts through
# --------------------------------------------------------------------------

def _ce(z0, z_ref: float) -> np.ndarray:
    """Neutral bulk exchange coefficient. `build_surface_roughness.py`'s own."""
    return KARMAN ** 2 / np.log(z_ref / np.maximum(np.asarray(z0), 1e-6)) ** 2


def _derived_z0(z0_cover, slope_variance, spacing_m):
    """Total roughness from the builder's derived scheme, on given cover values.

    The orographic term is a CELL statistic in both arms of the gap -- subgrid
    slope has no per-region value at the support this measures -- and the cover
    roughness is what the two arms disagree about. So the same slope variance
    and the same spacing are handed to both, and the difference between them is
    the reduction order and nothing else.

    Every function called here is `build_surface_roughness.py`'s own. The
    bisection this replaced was the builder's solve reproduced, which is what
    `failure-modes.md` class 17 forbids and what the derived scheme removed the
    need for: there is nothing left to solve.
    """
    lam = RESOLVED_BAND_MULTIPLE * np.asarray(spacing_m, dtype=np.float64)
    z0_cover = np.maximum(np.asarray(z0_cover, dtype=np.float64), 1e-6)
    h_m = pressure_scale_height(lam, z0_cover)
    l_in = inner_layer_depth(lam, z0_cover)
    z_m = np.maximum(h_m, lam * np.sqrt(2.0 * slope_variance) / np.pi)
    ca, cmd, _ = orographic_form_drag(slope_variance, z_m, h_m, l_in, z0_cover)
    return np.maximum(effective_length(ca + cmd, z_m), z0_cover)


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

    # The subgrid slope, and the mesh spacing it belongs to: the two inputs the
    # derived orographic term takes. Both are cell statistics, both are handed
    # unchanged to each arm.
    grad2 = np.tan(np.radians(mesh.local_slope_deg.astype(np.float64))) ** 2
    t2 = np.bincount(cell[land], weights=area[land] * 0.5 * grad2[land],
                     minlength=ncell)
    theta2 = np.where(have, t2 / np.maximum(w, 1e-30), 0.0)
    spacing_region = 1000.0 * np.sqrt(2.0 / np.sqrt(3.0)) * np.sqrt(area)
    s1 = np.bincount(cell[land], weights=area[land] * spacing_region[land],
                     minlength=ncell)
    spacing = np.where(have, s1 / np.maximum(w, 1e-30), 1.0e4)

    z0_cell = _derived_z0(np.where(have, z0_surf_cell, DEFAULT_BARE_Z0_M),
                          theta2, spacing)
    z0_cell = np.where(have, z0_cell, 0.0)
    z0_region = _derived_z0(z0_surface[land], theta2[cell[land]],
                            spacing[cell[land]])

    # The instrument, computed FIRST because every arm below is judged against
    # it: the step's own bracket on `ce` at the land-mean roughness, which is
    # the range `build_surface_roughness.py` already reports because `z_ref` is
    # not known before a climatology exists.
    ce_bracket = [float(_ce(np.array([np.average(z0_cell[have], weights=w[have])]),
                            reference_height_m(t, config))[0])
                  for t in CE_BRACKET_K]
    instrument = abs(ce_bracket[1] - ce_bracket[0]) / min(ce_bracket)

    out = {"derived_land_mean_z0_m": round(
               float(np.average(z0_cell[have], weights=w[have])), 6),
           "cover_land_mean_z0_m": round(
               float(np.average(z0_surf_cell[have], weights=w[have])), 6),
           "slope_variance_land_median": float(np.median(theta2[have])),
           "mesh_spacing_m_land_mean": round(
               float(np.average(spacing[have], weights=w[have])), 1),
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
    """How much of a cross-rung difference in roughness is calibration.

    It used to be most of it. `build_surface_roughness.py` solved a free
    orographic coefficient at every resolution so the land mean landed on a
    reference, so two rungs differed by their terrain and by their calibration
    at once and finding 8 of the audit could not be answered from them.

    The derived scheme has no coefficient to solve. Its orographic term is a
    property of the land and of the mesh the land is sampled on, and the grid
    rung enters only through which mesh regions fall in which cell, so the
    answer is expected to be flat across the ladder rather than merely
    comparable. That is the claim this arm now measures, and it can fail: a
    derived land mean that moved with the rung would say the scheme is still
    carrying the support inside it.
    """
    ref = per_rung[reference]["derived_land_mean_z0_m"]
    rows = {}
    for rung, arm in per_rung.items():
        rows[rung] = {
            "aggregate_then_process_land_mean_z0_m": arm["derived_land_mean_z0_m"],
            "ratio_to_reference": round(arm["derived_land_mean_z0_m"] / ref, 6)
                                  if ref else None,
            "cover_land_mean_z0_m": arm["cover_land_mean_z0_m"],
            "slope_variance_land_median": arm["slope_variance_land_median"],
            "subgrid_stdev_m_land_median": arm["subgrid_stdev_m_land_median"],
        }
    spread = [r["ratio_to_reference"] for r in rows.values()
              if r["ratio_to_reference"] is not None]
    return {"reference_rung": reference,
            "namelist_fallback_m": EXOPLASIM_DZ0LAND_M,
            "note": "the derived orographic term carries no solved coefficient, "
                    "so a difference between two rungs is terrain and "
                    "aggregation alone. The land means here are this arm's "
                    "AGGREGATE-THEN-PROCESS side -- the cover roughness is a "
                    "linear mean of lengths, which is the reduction the gap "
                    "exists to measure against -- so they are not the field's "
                    "own land mean, which build_surface_roughness.py reduces in "
                    "ce at the blending height and reports per rung. "
                    "`namelist_fallback_m` is the model's uniform default over "
                    "land and is neither a target nor an anchor.",
            "ratio_spread": round(max(spread) - min(spread), 6) if spread else None,
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

        depth = maximum_depth * P / (P + E),  E ~ erodibility

    which is convex in the erodibility, so the mean of the depths is above the
    depth at the mean. Erodibility is a per-region property of the rock the
    region is made of, so both orders are well defined and the comparison is
    the operator's and nothing else's.

    Everything the law needs except the erodibility collapses into ONE ratio,

        rho = E(cell mean erodibility) / P

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

def _critical_lapse(land_mean, t0: float, ceiling: float) -> float | None:
    """The lapse rate at which the land-mean gap reaches the bar, K per km.

    THIS IS NOT A NEW BAR. The bar is unchanged; this reports where it sits on
    the one axis the whole verdict turns on, so that a climatology settles the
    arm with a single comparison instead of a re-measurement. The gap rises
    monotonically with the lapse rate from zero at an isothermal column, so a
    bisection is exact to the tolerance stated.

    `None` means the bar is not reached anywhere this planet can lapse: the
    reduction is admissible at that rung whatever the climatology turns out to
    be, which is the one outcome that is a PASS rather than a deferral.
    """
    if land_mean(ceiling, t0) <= SATURATION_BAR:
        return None
    lo, hi = 0.0, ceiling
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if land_mean(mid, t0) > SATURATION_BAR:
            hi = mid
        else:
            lo = mid
    return round(0.5 * (lo + hi), 4)


def saturation_arm(mesh: Export, cell, ncell, land, area, config) -> dict:
    """What a mean-elevation cell costs in saturation vapour pressure.

    `build_boundary_conditions.py` averages elevation over a cell's land and
    hands the model one height. Temperature is LINEAR in height through a lapse
    rate, so that reduction is exact for temperature and the Jensen term is
    zero. Saturation vapour pressure is not: the model's Magnus-Teten function
    is convex in temperature over the whole liquid-water span, so the cell's own
    spread of heights raises the mean saturation deficit above the one computed
    at the mean height. That is the term a single elevation per cell discards,
    and evaporation is what consumes it.

    The function evaluated is the MODEL'S, phase switch included, because the
    quantity being measured is a curvature and the two functions do not have
    the same one. At the cold corner a cell's high ground falls below `tmelt`,
    where the model's coefficients are steeper than the liquid ones.

    Bracketed twice over, and both brackets are reported: the lapse rate,
    because which one applies is a property of a climatology this step does not
    have, and the reference air temperature, over the span in which surface
    water is liquid. The lapse ceiling is this planet's dry adiabat and the
    floor is declared with nothing under it, both per the criterion above.
    """
    coeff = magnus_teten()
    bracket = lapse_bracket_k_per_km(config)
    elev_km = mesh.field("elevation_km").astype(np.float64)
    a = area[land]
    c = cell[land]
    w = np.bincount(c, weights=a, minlength=ncell)
    have = w > 0
    m1 = np.bincount(c, weights=a * elev_km[land], minlength=ncell)
    mean = np.where(have, m1 / np.maximum(w, 1e-30), 0.0)
    anomaly = elev_km[land] - mean[c]
    land_area = w[have]

    def gap_field(lapse: float, t0: float) -> np.ndarray:
        # e_sat at each region's own height, against e_sat at the cell mean
        # height. `t0` is the temperature the mean height sits at.
        ratio = model_e_sat(t0 - lapse * anomaly, coeff) / model_e_sat(t0, coeff)
        num = np.bincount(c, weights=a * ratio, minlength=ncell)
        return np.where(have, num / np.maximum(w, 1e-30), 1.0) - 1.0

    def land_mean(lapse: float, t0: float) -> float:
        g = gap_field(lapse, t0)[have]
        return float((g * land_area).sum() / land_area.sum())

    rows = {}
    for lapse in bracket:
        for t0 in REFERENCE_AIR_BRACKET_K:
            gap = gap_field(lapse, t0)
            rows[f"lapse{lapse:g}_t{t0}"] = {
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
            "lapse_bracket_k_per_km": [round(v, 4) for v in bracket],
            "lapse_ceiling_is": "this planet's dry adiabat, lib/lapse.py",
            "lapse_floor_is": "declared; nothing at this step bounds it",
            "reference_air_bracket_k": list(REFERENCE_AIR_BRACKET_K),
            "critical_lapse_k_per_km": {
                str(t0): _critical_lapse(land_mean, t0, bracket[1])
                for t0 in REFERENCE_AIR_BRACKET_K},
            "saturation_function": "plasimmod.f90 Magnus-Teten, phase switch "
                                   "at tmelt, coefficients read from "
                                   "p_earth.f90",
            "magnus_teten_coefficients": coeff,
            "by_bracket_corner": rows}


# --------------------------------------------------------------------------
# arm 5: the support a length-carrying statistic is measured on
# --------------------------------------------------------------------------

def support_arm(config, rung: str, pedo) -> dict:
    """The within-cell elevation spread on this project's two builds, one grid.

    Both builds are the same planet at the same seed, so the ratio between them
    is what changing the region count alone does. The arm reports the spread
    over BOTH runs: the mesh spacing, which is what the statistic used to carry
    and which moves with the region count, and `catena.gradient_baseline_km`,
    the declared run `pedology/scripts/build_soil.py:subgrid_slope` now divides
    by. Which of the three shift mechanisms this one is, and the tests that
    separate them, are `analysis/subgrid_slope_support.py`; this arm only sizes
    the spread and what it is worth in metres of regolith through the catena
    divisor, which is the consuming step.
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
        spread_m = np.sqrt(np.maximum(v, 0.0))
        tanb = spread_m / spacing_m
        baseline_m = float(pedo["catena"]["gradient_baseline_km"]) * 1000.0
        tanb_declared = spread_m / baseline_m
        have = n > 1
        out[name] = {
            "regions": int(mesh.n_regions),
            "mesh_spacing_km": round(spacing_m / 1000.0, 4),
            "land_median_tan_beta": round(float(np.median(tanb[have])), 6),
            "land_median_tan_beta_declared_run": round(
                float(np.median(tanb_declared[have])), 6),
            "land_median_stdev_m": round(float(np.median(np.sqrt(
                np.maximum(v[have], 0.0)))), 3),
            "land_median_catena_divisor": round(
                float(np.median(1.0 / (1.0 + s * tanb[have]))), 6),
        }
    if all("absent" not in v for v in out.values()):
        a, b = out["precarve-craton"], out["precarve-craton-10m"]
        out["ratio_10m_over_2p5m"] = {
            "tan_beta": round(b["land_median_tan_beta"] / a["land_median_tan_beta"], 4),
            "tan_beta_declared_run": round(
                b["land_median_tan_beta_declared_run"]
                / a["land_median_tan_beta_declared_run"], 4),
            "stdev_m": round(b["land_median_stdev_m"] / a["land_median_stdev_m"], 4),
            "catena_divisor": round(
                b["land_median_catena_divisor"] / a["land_median_catena_divisor"], 4),
        }
    return {"rung": rung, "slope_transport": s, "builds": out}


# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# arm 7: the staged dry and saturated albedo pair through the model's mixing
# --------------------------------------------------------------------------

# The dispositions world-yxor names, priced against the SAME truth so that what
# separates them is the staging and nothing else. `as_built` is the defect,
# `staged_in_transform` is the repair the audit pre-registered and refuted, and
# `three_point_knee` is the one-extra-field form of the model change. A
# candidate is VIABLE when its absorbed-shortwave gap is inside the storage
# tolerance at every saturation the modelled column reaches; a candidate that is
# not is refuted, and that is a result rather than a reason to try another.
CANDIDATE_REPAIRS = ("as_built", "staged_in_transform", "three_point_knee")

LANDMOD = PROJECT_ROOT / "vendor/exoplasim/exoplasim/plasim/src/landmod.f90"
_DRHSFULL = re.compile(r"^\s*real\s*::\s*drhsfull\s*=\s*([0-9.eEdD+-]+)", re.M)


def evaporation_knee_fill_fraction(path: Path = LANDMOD) -> float:
    """`drhsfull`, READ from the model source rather than copied beside it.

    It is the fill fraction of the surface layer above which the evaporation
    limiter's wetness factor reaches one, and it is landmod's own: nothing routes
    it through a namelist, so a value written into this file would be a
    declaration nothing re-derives and nothing objects to when the two disagree.

    It matters here because it is where the modelled soil's albedo is evaluated
    most often rather than where a gap happens to be largest. The limiter's knee
    is the transition the limiter exists to describe, so a soil that wets and
    dries through the limiter's ceiling crosses it every cycle.
    """
    text = path.read_text(encoding="utf-8")
    match = _DRHSFULL.search(text)
    if match is None:
        raise SystemExit(
            f"{path} declares no `real :: drhsfull = ...`; the evaporation "
            "knee cannot be read and must not be guessed")
    return float(match.group(1).replace("d", "e").replace("D", "e"))


def _km_inverse(r):
    """The Kubelka-Munk inverse `R = 1 + r - sqrt(r^2 + 2r)`.

    `sadeghi_mix` applies it at the end of the mixing; here it is needed on its
    own to build the staged-in-the-transform candidate, which is what the note's
    refuted repair actually stages.
    """
    r = np.asarray(r, dtype=np.float64)
    return 1.0 + r - np.sqrt(r * r + 2.0 * r)


def _three_point_mix(a_dry, a_knee, a_wet, s, s_knee, sigma):
    """Sadeghi's mixing in two segments through a third staged field.

    THE FORM, and it is a declared model form rather than a fit. Each of the
    three staged fields is an exact area mean of a per-region quantity -- the dry
    albedo, the saturated albedo, and the cell's own mixed albedo at the knee --
    so the composition is exact at all three saturations by construction, where
    the two-field form is exact at two. Between them the model runs the same
    Sadeghi curve on the saturation rescaled inside the segment it is in.

    Nothing here is chosen to make a comparison come out. The only choice is
    WHERE the third point sits, and that is the evaporation limiter's own knee,
    read from `landmod.f90` and mapped through the declared saturation
    endpoints.
    """
    s = np.asarray(s, dtype=np.float64)
    lower = s <= s_knee
    sub = np.where(lower, s / s_knee, (s - s_knee) / (1.0 - s_knee))
    return np.where(lower,
                    sadeghi_mix(a_dry, a_knee, sub, sigma),
                    sadeghi_mix(a_knee, a_wet, sub, sigma))


def _region_albedo_pair(mesh: Export, config: dict):
    """Per-region dry and saturated albedo, per band, as the builder forms them.

    Every expression here is `build_surface_albedo.py`'s own, imported rather
    than restated: the band shapes, the wetting ratios and the lithology
    overrides all come from the step that owns them. What this does NOT do is
    the reduction, which is the whole point of the arm; both orders of it are
    taken by the caller.

    The vegetation, lake and evaporite repaints are deliberately absent. Each is
    applied per region BEFORE the reduction in the builder, so they change which
    class a region carries and therefore the size of the within-cell spread, and
    they do not change the ORDER of the operation. Leaving them out measures the
    reduction on the lithology alone, which is the population that exists on
    terrain with no climate run and no biosphere.
    """
    rock = mesh.field("substrate_class").astype(int)
    region_albedo = mesh.rock_albedo.astype(np.float64).copy()

    overrides = (config.get("model") or {}).get("lithology_albedo_overrides") or {}
    applied = {}
    for code, spec in overrides.items():
        value = float(spec["albedo"] if isinstance(spec, dict) else spec)
        rid = _rock_id(mesh.root, code)
        sel = rock == rid
        if sel.any():
            applied[code] = {"rock_id": rid, "albedo": value,
                             "regions": int(sel.sum())}
            region_albedo[sel] = value

    z1, z2 = (float(v) for v in band_fractions())
    rho = np.empty(region_albedo.shape, dtype=np.float64)
    for rid, value in rock_band_ratios(ROCK_BANDS, mesh.root).items():
        rho[rock == rid] = value
    shape1, shape2 = band_shapes(rho, z1, z2)

    rock_wet, wet_report = rock_wetting_ratios(WETTING, mesh.root, WETTING_ARM)
    wet1 = np.empty(region_albedo.shape, dtype=np.float64)
    wet2 = np.empty(region_albedo.shape, dtype=np.float64)
    for rid, (w1, w2) in rock_wet.items():
        wet1[rock == rid] = w1
        wet2[rock == rid] = w2

    dry = (region_albedo * shape1, region_albedo * shape2)
    sat = (dry[0] * wet1, dry[1] * wet2)
    return rock, dry, sat, {"band_fractions": [z1, z2],
                            "lithology_albedo_overrides": applied,
                            "wetting": wet_report}


def _one_class_cells(cell, ncell: int, rock, land) -> np.ndarray:
    """True where every land region in the cell carries the same substrate class."""
    order = np.argsort(cell[land], kind="stable")
    c = cell[land][order]
    r = rock[land][order]
    starts = np.flatnonzero(np.r_[True, c[1:] != c[:-1]])
    ends = np.r_[starts[1:], c.size]
    mixed = np.zeros(ncell, dtype=bool)
    for lo, hi in zip(starts, ends):
        block = r[lo:hi]
        mixed[c[lo]] = bool((block != block[0]).any())
    return ~mixed


def band_downward_shortwave(grid: tuple[int, int]) -> dict | None:
    """Annual-mean downward surface shortwave per band, on the accepted baseline.

    Returns None where the accepted climatology is not on this grid. There is no
    fallback onto another rung: a plausible flux from a grid the albedo field was
    not built on is worse than no flux, which is the disposition
    `analysis/coastline_flux_bracket.py` already takes for the same reason.

    The two axes are joined by INDEX and never by label. `require_same_rows` is
    the door for an axis read off a CLIMATOLOGY: netCDF stores it as float32, so
    the constructed Gauss-Legendre nodes and the axis on disk agree to a few
    parts in a million and no closer, and `require_gaussian_rows`'s float64 bar
    refuses a correct axis for that reason alone. The columns are the identity,
    because the export and the model label the same columns differently and only
    the labels differ. `CLAUDE.md` rule 3.
    """
    from netCDF4 import Dataset
    clim = best_available_climatology()
    if clim is None or not Path(clim.path).is_file():
        return None
    with Dataset(clim.path) as ds:
        shape = (len(ds.dimensions["lat"]), len(ds.dimensions["lon"]))
        if shape != grid:
            return None
        require_same_rows(gaussian_latitudes(grid[0]),
                          np.asarray(ds["lat"][:], dtype=np.float64),
                          "the accepted climatology against this rung")
        centres = np.asarray(ds["time"][:], dtype=np.float64)
        bands = [annual_mean(np.asarray(ds[n][:], dtype=np.float64), centres).ravel()
                 for n in ("rsds1", "rsds2")]
    return {"stage": clim.stage, "path": str(clim.path), "bands": bands}


def wet_albedo_arm(mesh: Export, cell, ncell, land, area, config,
                   grid: tuple[int, int]) -> dict:
    """The staged pair reduced then mixed, against mixed then reduced.

    `build_surface_albedo.py` stages a DRY band pair and a SATURATED band pair,
    each an area-weighted mean over the cell's own lithology, and
    `landmod.f90:wetalb` mixes them per cell through `wet_soil_albedo`, which is
    Sadeghi, Jones and Philpot (2015): linear in the Kubelka-Munk transform
    `r = (1-R)^2 / (2R)` and therefore NOT in the albedo. So the model computes
    the mixing of the cell's mean rock, where what the cell owes the atmosphere
    is the area mean of the mixings of its own rocks.

    The saturation axis is not swept over a guess. Its two endpoints are
    `soil_albedo_moisture.saturation_at_empty_layer` and
    `saturation_at_full_layer`, which `run_exoplasim.py` already checks against
    the compiled defaults, so the sweep covers exactly the range the modelled
    land column reaches and no more.

    The bar is the accepted baseline's 0.12 W m-2 state-storage tolerance. An
    albedo error reaches the surface energy balance as `delta_alpha * S_down`,
    so the arm weights the per-cell gap by the accepted baseline's own downward
    shortwave in each band and by the cell's land area and divides by the area
    of the planet. That number is available on the rung the accepted
    climatology carries and on no other, and the arm says which.
    """
    moisture = config["surface"]["soil_albedo_moisture"]
    s_lo = float(moisture["saturation_at_empty_layer"])
    s_hi = float(moisture["saturation_at_full_layer"])
    sigma = (float(moisture["shape_sigma_band1"]),
             float(moisture["shape_sigma_band2"]))
    sweep = tuple(round(s_lo + (s_hi - s_lo) * k / 8.0, 6) for k in range(9))

    rock, dry, sat, meta = _region_albedo_pair(mesh, config)
    one_class = _one_class_cells(cell, ncell, rock, land)
    shortwave = band_downward_shortwave(grid)

    # THE THIRD POINT'S SATURATION, derived rather than chosen. `drhsfull` is
    # read from `landmod.f90` and mapped onto saturation by the same two declared
    # endpoints `wetalb` maps a fill fraction through, so the knee is the image
    # of a model constant under a config declaration and neither is written here.
    fill_knee = evaporation_knee_fill_fraction()
    s_knee = s_lo + fill_knee * (s_hi - s_lo)
    if not s_lo < s_knee < s_hi:
        raise SystemExit(
            f"the evaporation knee maps to a saturation of {s_knee}, outside "
            f"the declared range {s_lo} to {s_hi}; a third staged point outside "
            "the range the column reaches is not a repair")

    land_area = np.bincount(cell[land], weights=area[land], minlength=ncell)
    covered = land_area > 0
    total_area = np.bincount(cell, weights=area, minlength=ncell)
    planet_area = float(total_area.sum())

    worst_pure = 0.0
    bands: dict[str, dict] = {}
    watt_gap = {s: np.zeros(ncell) for s in sweep}
    # Every candidate form is priced in the SAME terms as the defect and on the
    # same truth, so what separates them is the staging and nothing else.
    candidate_watts = {name: {s: np.zeros(ncell) for s in sweep}
                       for name in CANDIDATE_REPAIRS}
    candidate_albedo: dict[str, dict[str, list]] = {
        name: {} for name in CANDIDATE_REPAIRS}
    for b, (d_r, w_r, sig) in enumerate(zip(dry, sat, sigma), start=1):
        a_dry, _, cov = cell_mean(cell, ncell, area, d_r, land)
        a_wet, _, _ = cell_mean(cell, ncell, area, w_r, land)
        # The two extra stagings the candidates need, each a reduction of a
        # per-region quantity and neither a free parameter.
        r_dry_bar, _, _ = cell_mean(cell, ncell, area,
                                    kubelka_munk_transform(d_r), land)
        r_wet_bar, _, _ = cell_mean(cell, ncell, area,
                                    kubelka_munk_transform(w_r), land)
        a_dry_t, a_wet_t = _km_inverse(r_dry_bar), _km_inverse(r_wet_bar)
        a_knee, _, _ = cell_expectation(
            cell, ncell, area,
            lambda values, w=w_r, sig=sig: sadeghi_mix(values, w, s_knee, sig),
            d_r, land)
        candidate_rows: dict[str, list] = {name: [] for name in CANDIDATE_REPAIRS}
        rows, stack = [], []
        for s in sweep:
            aggregate = sadeghi_mix(a_dry, a_wet, s, sig)
            # `cell_expectation` applies the law to the WHOLE per-region array
            # and averages afterwards, so the saturated partner closes over the
            # same region order and the alignment is the operator's own.
            expectation, _, _ = cell_expectation(
                cell, ncell, area,
                lambda values, w=w_r, s=s, sig=sig: sadeghi_mix(values, w, s, sig),
                d_r, land)
            gap = np.where(cov, expectation - aggregate, 0.0)
            stack.append(gap)
            for name, form in (
                    ("as_built", aggregate),
                    ("staged_in_transform", sadeghi_mix(a_dry_t, a_wet_t, s, sig)),
                    ("three_point_knee", _three_point_mix(
                        a_dry, a_knee, a_wet, s, s_knee, sig))):
                cgap = np.where(cov, expectation - form, 0.0)
                if shortwave is not None:
                    candidate_watts[name][s] = (candidate_watts[name][s]
                                                + cgap * shortwave["bands"][b - 1])
                cw, cg = land_area[cov], cgap[cov]
                candidate_rows[name].append({
                    "saturation": s,
                    "land_mean_gap": float((cg * cw).sum() / cw.sum()),
                    "maximum_absolute_gap": float(np.abs(cg).max()),
                })
            pure = cov & one_class
            if pure.any():
                worst_pure = max(worst_pure, float(np.abs(gap[pure]).max()))
            if shortwave is not None:
                watt_gap[s] = watt_gap[s] + gap * shortwave["bands"][b - 1]
            w = land_area[cov]
            g = gap[cov]
            rows.append({
                "saturation": s,
                "land_mean_gap": float((g * w).sum() / w.sum()),
                "land_mean_aggregate_then_process": float(
                    (aggregate[cov] * w).sum() / w.sum()),
                "per_cell_gap": {str(p): float(np.percentile(g, p))
                                 for p in (1, 50, 99)},
                "maximum_absolute_gap": float(np.abs(g).max()),
            })
        for name in CANDIDATE_REPAIRS:
            candidate_albedo[name][f"band{b}"] = candidate_rows[name]
        # The pre-registered sign invariant, reported as what can falsify it:
        # the number of MIXED cells whose gap takes both signs across the
        # interior of the saturation sweep.
        interior = np.stack([g for g, s in zip(stack, sweep) if s_lo < s < s_hi])
        both = (interior > 0).any(axis=0) & (interior < 0).any(axis=0)
        bands[f"band{b}"] = {
            "shape_sigma": sig,
            "by_saturation": rows,
            "cells_whose_gap_changes_sign_along_saturation": int(
                (both & cov & ~one_class).sum()),
        }

    candidates: dict[str, dict] = {}
    for name in CANDIDATE_REPAIRS:
        row: dict = {"albedo_gap_by_band": candidate_albedo[name]}
        if shortwave is not None:
            watts = {str(s): float((candidate_watts[name][s] * land_area).sum()
                                   / planet_area) for s in sweep}
            peak = max(abs(v) for v in watts.values())
            row["absorbed_shortwave_w_m2"] = watts
            row["peak_absolute_w_m2"] = peak
            row["saturations_past_the_bar"] = [
                s for s in watts if abs(watts[s]) > STORAGE_TOLERANCE_W_M2]
            row["inside_the_bar_everywhere"] = peak <= STORAGE_TOLERANCE_W_M2
        candidates[name] = row

    result = {
        "saturation_endpoints": [s_lo, s_hi],
        "saturation_sweep": list(sweep),
        "wetting_arm": WETTING_ARM,
        "candidate_repairs": {
            "why": ("world-yxor offers a model change or a declared bracket, and "
                    "the choice cannot be taken until the model change is priced: "
                    "the repair this audit pre-registered for the same defect was "
                    "measured wrong by two orders, so a proposed form is measured "
                    "before it is adopted"),
            "criterion": ("VIABLE when the absorbed-shortwave gap is inside "
                          f"{STORAGE_TOLERANCE_W_M2} W m-2 at every saturation in "
                          "the declared range; the criterion is the same bar the "
                          "defect is measured against and predates the result"),
            "evaporation_knee": {
                "fill_fraction": fill_knee,
                "source": "landmod.f90, real :: drhsfull",
                "saturation": s_knee,
                "mapped_by": ("config/planet.yaml soil_albedo_moisture "
                              "saturation_at_empty_layer/saturation_at_full_layer, "
                              "which is the map wetalb applies"),
            },
            "forms": candidates,
        },
        "cells_with_land": int(covered.sum()),
        "cells_of_one_substrate_class": int((covered & one_class).sum()),
        "one_class_control": {
            "worst_absolute_gap": worst_pure,
            "tolerance": ONE_CLASS_TOLERANCE,
            "passes": worst_pure <= ONE_CLASS_TOLERANCE,
        },
        "bands": bands,
        "inputs": meta,
    }
    if shortwave is None:
        result["absorbed_shortwave"] = {
            "available": False,
            "reason": "the accepted climatology is not on this rung, and there "
                      "is no fallback onto another one",
        }
    else:
        result["absorbed_shortwave"] = {
            "available": True,
            "stage": shortwave["stage"],
            "bar_w_m2": STORAGE_TOLERANCE_W_M2,
            "bar_source": "config/partial_surface.yaml, the accepted baseline "
                          "state-storage tolerance the tile operator was "
                          "selected against",
            "global_mean_by_saturation": {
                str(s): float((watt_gap[s] * land_area).sum() / planet_area)
                for s in sweep},
            "land_mean_downward_shortwave_w_m2": [
                float((band * land_area).sum() / land_area.sum())
                for band in shortwave["bands"]],
        }
    return result


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path,
                    default=PROJECT_ROOT / "config" / "planet.yaml")
    ap.add_argument("--rung", action="append", default=None,
                    help="ladder rung to measure; repeatable, defaults to "
                         "every rung the configured build carries an export for")
    ap.add_argument("--reference-rung", default="T42",
                    help="rung the roughness recalibration is compared against")
    ap.add_argument("--build", default=None,
                    help="build to measure; the configured one by default. "
                         "Every arm here is terrain, so a build with a payload "
                         "is measurable whether or not it is the active one")
    ap.add_argument("--output", type=Path, default=OUTPUT)
    args = ap.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    pedo = yaml.safe_load(
        (PROJECT_ROOT / "pedology" / "config" / "pedogenesis.yaml")
        .read_text(encoding="utf-8"))

    build_root = (builds.build_root(config) if args.build is None
                  else PROJECT_ROOT / "source" / args.build)
    mesh = Export(builds.mesh_export_at(build_root))
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
            "saturation_deficit": saturation_arm(mesh, cell, ncell, land,
                                                 area, config),
            "wet_albedo_mixing": wet_albedo_arm(mesh, cell, ncell, land, area,
                                                config, (nlat, nlon)),
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
            "saturation_deficit": f"{SATURATION_BAR} of e_sat, declared; "
                                  "unchanged since the first measurement "
                                  "crossed it",
            "wet_albedo_mixing": f"{STORAGE_TOLERANCE_W_M2} W m-2, the accepted "
                                 "baseline state-storage tolerance the tile "
                                 "operator was selected against, reached by an "
                                 "albedo error as delta_alpha * S_down",
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
