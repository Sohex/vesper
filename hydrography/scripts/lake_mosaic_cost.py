"""What the hydrologic mosaic costs at the climate-grid crossing, and what the
seasonal lake-energy omission is worth once a baseline climatology exists.

Worldbuilding. Vesper is an invented super-Earth. Everything here is about the
simulation of its land surface: the modelled lake set solved by
`surface_water.py`, the bucket capacity `landmod.f90` runs on, and the baseline
climatology's own surface fields.

Two rows, one instrument, because both read the same three artifacts.

THE CAPACITY THE FLUX ARMS STAND ON IS THE ONE THE CLIMATOLOGY'S OWN RUN
INTEGRATED, not the one staged now. The wetness factor is a function of the
soil water and the capacity together, so the two have to come from one state:
`exoplasim/inputs/<rung>/` holds the field the NEXT run will read, and a
restage is an ordinary event. On this build the two differ on 1598 of 2048
cells, and pairing the run's soil water with the staged bucket puts the soil
water above the bucket on 55 land cells, which a clipped store cannot do.
`lib/provenance.py:run_surface_field` is the door onto the field a run
consumed; `staged_surface_field` is the door onto the field the next run will
read, and the perturbation arm keeps that one because it is about the next
staging. world-q3p9.

HYD-21 asks for the seasonal lake-energy omission to be BOUNDED before a lake
implementation is chosen. `notes/audits/lake-energy-omission-bound.md` bounds
it in closed form and leaves one term open, because it was written while
`config/planet.yaml` carried `baseline_climatology: null`: freezing dominates
the other two terms, and whether a lake freezes is a climatology question. A
baseline climatology now exists, so this closes it.

LSHY-6 asks the native-mesh to climate-grid crossing to keep the upland-soil
and lake/playa tiles distinct. `build_spatial_support.py` does keep them
distinct; the loss is at the CONSUMER, where the two-tile capacity distribution
is collapsed to one scalar before two convex laws read it. This measures both
the collapse that is happening today and the one that would remain if the
blend were staged.

THE INSTRUMENT IS BOUNDED RATHER THAN DIVIDED, and the reason is in the
numbers. The wetness factor is `min(1, w / (0.4 C))` at `landcolumn.f90:370`
and the latent heat flux is that factor times a potential flux. Recovering the
potential flux as `hfls / beta` looks exact and is not usable: `beta` is under
0.05 on nearly half of this planet's land bins, and the quotient reaches 1e17
W m-2. The potential flux is therefore CAPPED at the surface energy actually
available, `rss + rls`, and the result is reported as a bracket whose upper end
assumes every cell could put its whole net radiation into evaporation. A bound
that holds at the upper end holds.

Reads only. Nothing here regenerates a tracked artifact.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import _paths  # noqa: F401
from _paths import ANALYSIS, DATA, PROJECT_ROOT

import netCDF4 as nc  # noqa: E402
import yaml  # noqa: E402
import provenance  # noqa: E402
from paths import climatology_path, rel  # noqa: E402

# From notes/audits/lake-energy-omission-bound.md, both derived there in closed
# form from the model's own constants and this planet's orbital period.
SURROGATE_HEAT_CAPACITY_J_M2_K = 3.295e6
ICE_LATENT_HEAT_J_M2_PER_M = 3.06e8

FRESH_FREEZING_K = 273.15
# The deepest a freezing point goes for a closed-basin chloride brine: the
# NaCl eutectic. This world's lake chemistry is not declared, so the depression
# is used as a BRACKET WIDTH and never as this planet's freezing point.
EUTECTIC_DEPRESSION_K = 21.1

# landcolumn.f90:370. The wetness factor saturates at this share of capacity.
DRHSFULL = 0.4

# config/partial_surface.yaml's materiality test: the accepted baseline run's
# own state-storage tolerance, which that file already selected the tile
# operator against and which the land audit's albedo arm reports in.
BAR_W_M2 = 0.12


def read_sra_field(path: Path) -> np.ndarray:
    """The values of a one-code `.sra`, header dropped."""
    values: list[float] = []
    for line in [ln for ln in path.read_text().split("\n") if ln.strip()][1:]:
        values.extend(float(x) for x in line.split())
    return np.asarray(values, float)


def wetness(water: np.ndarray, capacity: np.ndarray) -> np.ndarray:
    """`drhs`, landcolumn.f90:370, on a capacity that may be zero."""
    safe = np.where(capacity > 0.0, capacity, np.nan)
    return np.minimum(1.0, water / (DRHSFULL * safe))


def area_weighted_quantile(values, weights, probabilities):
    order = np.argsort(values)
    v, w = np.asarray(values)[order], np.asarray(weights)[order]
    if w.sum() <= 0:
        return {p: float("nan") for p in probabilities}
    cumulative = np.cumsum(w) / w.sum()
    return {p: float(np.interp(p, cumulative, v)) for p in probabilities}


def share(weights: np.ndarray, mask: np.ndarray) -> float:
    total = weights.sum()
    return float(100.0 * weights[mask].sum() / total) if total > 0 else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", default=None,
                    help="build name; defaults to config's source_build")
    ap.add_argument("--rung", default="T21",
                    help="ladder rung; only rungs with a support file can be read")
    ap.add_argument("--out", type=Path, default=ANALYSIS / "lake_mosaic_cost.json")
    ap.add_argument("--selftest", action="store_true",
                    help="run the instrument's own identities and exit; needs "
                         "no build, no climatology and no staged field")
    args = ap.parse_args()

    if args.selftest:
        selftest()
        return

    config = yaml.safe_load((PROJECT_ROOT / "config" / "planet.yaml").read_text())
    build = args.build or config["source_build"]
    lake_dwmax_m = float(config["model"]["lake_dwmax_m"])

    support_path = DATA / build / f"support_exoplasim-{args.rung}.nc"
    if not support_path.exists():
        raise SystemExit(
            f"{support_path} is absent. The mosaic shares come from SPAT-3's "
            "support file; run hydrography/scripts/build_spatial_support.py.")

    # THE BASELINE AND NOT THE BEST AVAILABLE. The freezing term is a threshold
    # on the settled climate STATE, and a bootstrap climatology is the
    # terrain-only-field answer to a question about the surface fields. A step
    # that fell back to the bootstrap here would return a plausible freezing
    # share for a world with no soil, no lakes and no vegetation in its albedo.
    # `climatology_path` refuses instead, and carries `require_configured_grid`
    # with it, which the hand-rolled read of the config key does not.
    clim_path = climatology_path(root=PROJECT_ROOT)

    support = nc.Dataset(support_path)
    clim = nc.Dataset(clim_path)

    # CLAUDE.md rule 3. The two label the same columns differently, so the
    # correspondence is BY INDEX and no longitude is ever compared. The only
    # thing checked is that the two grids ARE the same grid.
    for axis in ("lat", "lon"):
        if support.dimensions[axis].size != clim.dimensions[axis].size:
            raise SystemExit(
                f"the support is {support.dimensions['lat'].size}x"
                f"{support.dimensions['lon'].size} and the climatology is "
                f"{clim.dimensions['lat'].size}x{clim.dimensions['lon'].size}. "
                "These are different rungs; there is no mapping between them "
                "here and matching by coordinate value is refused.")
    shape = (support.dimensions["lat"].size, support.dimensions["lon"].size)

    # TWO CAPACITY FIELDS, AND THEY ANSWER DIFFERENT QUESTIONS. The flux arms
    # need a wetness factor, which is a function of the soil water AND the
    # capacity, so the two have to be one state: the climatology's soil water
    # is in equilibrium with the bucket THAT run integrated, and pairing it
    # with whatever is staged now is a read across two iterations. The staged
    # field answers the other question -- what the NEXT run will read -- and
    # the perturbation arm below is about exactly that, so it keeps it.
    run_id = getattr(clim, "vesper_run_id", None)
    if not run_id:
        raise SystemExit(
            f"{rel(clim_path)} carries no vesper_run_id, so which run's bucket "
            "its soil water is in equilibrium with is unrecorded and the flux "
            "arms have no self-consistent state to stand on.")
    consumed = provenance.run_surface_field(229, str(run_id))
    staged = provenance.staged_surface_field(229, config, for_build=build)
    staged_record = json.loads(
        (PROJECT_ROOT / staged["path"]).with_name(
            Path(staged["path"]).stem + "_provenance.json").read_text())

    capacity_flat = read_sra_field(PROJECT_ROOT / consumed["path"])
    staged_flat = read_sra_field(PROJECT_ROOT / staged["path"])
    for label, flat in (("the run's", capacity_flat), ("the staged", staged_flat)):
        if flat.size < shape[0] * shape[1]:
            raise SystemExit(f"{label} 0229 field holds {flat.size} values, "
                             f"fewer than the {shape[0] * shape[1]} the "
                             f"{args.rung} grid needs.")
    capacity = capacity_flat[:shape[0] * shape[1]].reshape(shape)
    staged_capacity = staged_flat[:shape[0] * shape[1]].reshape(shape)
    if consumed["source_build"] != build:
        raise SystemExit(
            f"{run_id} ran on {consumed['source_build']} and this is reading "
            f"{build}. Rule 5: a cross-build read is declared by naming the "
            "build, never defaulted.")
    # An inert blend is a property of a GENERATOR INVOCATION, and the record of
    # that invocation lives beside the staged field rather than travelling with
    # the copy in the run directory. So this is checked for the field that will
    # be staged and is unrecorded for the field that ran. world-hl06.
    blend_is_inert = staged_record.get("lakes") is None

    lake_km2 = np.asarray(support["solved_lake_area_km2"][:], float)
    land_km2 = np.asarray(support["land_area_km2"][:], float)
    cell_km2 = np.asarray(support["cell_area_km2"][:], float)
    f_lake = np.asarray(support["f_solved_lake"][:], float)
    covered = np.asarray(support["region_count"][:]) > 0
    is_land = land_km2 > 0
    has_lake = lake_km2 > 0

    ts = np.asarray(clim["ts"][:], float)
    mint = np.asarray(clim["mint"][:], float)
    soil_water = np.asarray(clim["mrso"][:], float)
    latent = np.abs(np.asarray(clim["hfls"][:], float))
    net_radiation = (np.asarray(clim["rss"][:], float)
                     + np.asarray(clim["rls"][:], float))
    if "bin_weight" in clim.variables:
        bin_weight = np.asarray(clim["bin_weight"][:], float)
    else:
        bin_weight = np.full(clim.dimensions["time"].size, 1.0)
    bin_weight = bin_weight / bin_weight.sum()

    report: dict = {
        "title": "What the hydrologic mosaic costs at the crossing, and what "
                 "the seasonal lake-energy omission is worth",
        "measured_utc": datetime.now(timezone.utc).isoformat(),
        "source_build": build,
        "rung": args.rung,
        "terrain_hash": getattr(support, "terrain_hash", None),
        "baseline_climatology": rel(clim_path),
        "baseline_run": getattr(clim, "vesper_run_id", None),
        "baseline_orbits": int(getattr(clim, "climatology_orbit_count", 0)) or None,
        "capacity_the_flux_arms_stand_on": consumed["path"],
        "capacity_sha256": consumed["sha256"],
        "staged_capacity": staged["path"],
        "staged_capacity_sha256": staged["sha256"],
        "staged_capacity_is_the_one_that_ran": consumed["matches_staged"],
        "lake_dwmax_m": lake_dwmax_m,
        "lake_blend_is_inert_in_the_staged_field": blend_is_inert,
        "bar_w_m2": BAR_W_M2,
        "bar_source": "config/partial_surface.yaml, the accepted baseline "
                      "run's state-storage tolerance",
    }

    # ---- controls -------------------------------------------------------
    surface_share = np.asarray(support["surface_share"][:], float)
    f_barren = np.asarray(support["f_barren"][:], float)
    f_nonbarren = np.asarray(support["f_nonbarren"][:], float)
    mesh_km2 = np.asarray(support["mesh_area_km2"][:], float)
    parts = (land_km2 + np.asarray(support["ocean_area_km2"][:], float)
             + np.asarray(support["inland_water_area_km2"][:], float))
    overshoot = soil_water.max(axis=0) - capacity
    staged_overshoot = soil_water.max(axis=0) - staged_capacity
    report["controls"] = {
        "surface_share_partition_residual": float(
            np.abs(surface_share.sum(axis=2)[covered] - 1.0).max()),
        "barren_partition_residual": float(
            np.abs(f_barren[is_land] + f_nonbarren[is_land] - 1.0).max()),
        "area_partition_relative_residual": float(np.abs(
            (parts - mesh_km2) / np.where(mesh_km2 > 0, mesh_km2, 1.0))[covered].max()),
        "land_cells": int(is_land.sum()),
        "land_cells_whose_soil_water_exceeds_the_capacity_that_ran": int(
            (overshoot[is_land] > 1e-6).sum()),
        "worst_overshoot_m": float(overshoot[is_land].max()),
        "max_fill_fraction": float(
            (soil_water.max(axis=0)[is_land]
             / np.where(capacity[is_land] > 0, capacity[is_land], np.nan)).max()),
        "land_cells_whose_soil_water_exceeds_the_STAGED_capacity": int(
            (staged_overshoot[is_land] > 1e-6).sum()),
        "worst_staged_overshoot_m": float(staged_overshoot[is_land].max()),
        "note": "the partition residuals are the float32 storage precision of "
                "the shares, not a reduction error; the areas are float64 and "
                "close exactly. THE OVERSHOOT COUNT IS A CONSERVATION TEST "
                "WITH A RIGHT ANSWER and the right answer is zero: "
                "bucket_step at landcolumn.f90:75-89 sets the store to "
                "min(capacity, store + flux dt) at every timestep and the "
                "layered path clips each layer at its share of the same "
                "capacity, so a bin mean of clipped values cannot exceed the "
                "clip. The maximum fill fraction is the same statement from "
                "the other side: the bucket touches its capacity and does not "
                "pass it. The STAGED count is reported beside it because it is "
                "not zero -- the staged field is a later iteration than the "
                "run and its bucket is a different bucket. world-q3p9.",
    }

    # A BUCKET CANNOT EXCEED ITS CAPACITY, so this is a refusal and not a
    # caveat. The tolerance is the `.sra` text format's own five decimals; the
    # measured headroom's minimum is 7e-6 m, so a real violation clears it by
    # orders of magnitude and the format's rounding does not reach it.
    violating = int((overshoot[is_land] > 1e-5).sum())
    if violating:
        raise SystemExit(
            f"the climatology's soil water exceeds the bucket capacity "
            f"{run_id} ran on, on {violating} of {int(is_land.sum())} land "
            f"cells, worst {float(overshoot[is_land].max()):.5f} m. "
            f"landcolumn.f90:75-89 clips the store at the capacity every "
            f"timestep, so this cannot happen for the run's own field and "
            f"means the pairing is wrong rather than the physics.")

    # ---- HYD-21: the freezing threshold ---------------------------------
    seasonal_min_warm = ts.min(axis=0)     # bin means understate the minimum
    seasonal_min_cold = mint.min(axis=0)   # within-bin minima overstate it
    seasonal_range = ts.max(axis=0) - ts.min(axis=0)
    equivalent_ice_m = (seasonal_range * SURROGATE_HEAT_CAPACITY_J_M2_K
                        / ICE_LATENT_HEAT_J_M2_PER_M)

    freezing = {}
    for name, weight in (("solved_lake_area", lake_km2), ("all_land_area", land_km2)):
        freezes = seasonal_min_warm < FRESH_FREEZING_K
        freezing[name] = {
            "share_below_the_fresh_freezing_point_bin_mean_end": share(weight, freezes),
            "share_below_the_fresh_freezing_point_within_bin_end": share(
                weight, seasonal_min_cold < FRESH_FREEZING_K),
            "share_below_the_eutectic": share(
                weight, seasonal_min_warm < FRESH_FREEZING_K - EUTECTIC_DEPRESSION_K),
            "share_inside_the_depression_band": share(
                weight, freezes & (seasonal_min_warm
                                   >= FRESH_FREEZING_K - EUTECTIC_DEPRESSION_K)),
        }
    population = has_lake
    freezing["equivalent_ice_thickness_m_over_lake_area"] = area_weighted_quantile(
        equivalent_ice_m[population], lake_km2[population], [0.05, 0.5, 0.95, 0.99])
    freezing["equivalent_ice_thickness_m_over_freezing_lake_area"] = area_weighted_quantile(
        equivalent_ice_m[population & (seasonal_min_warm < FRESH_FREEZING_K)],
        lake_km2[population & (seasonal_min_warm < FRESH_FREEZING_K)],
        [0.05, 0.5, 0.95])
    freezing["max_equivalent_ice_thickness_m"] = float(equivalent_ice_m[population].max())
    freezing["max_seasonal_range_K_over_lake_cells"] = float(seasonal_range[population].max())
    freezing["lake_area_share_whose_range_exceeds_the_93_K_figure"] = share(
        lake_km2, has_lake & (seasonal_range > 93.0))
    freezing["definition"] = (
        "equivalent ice thickness is the depth of lake ice whose latent heat "
        "of freezing equals the heat the modelled soil surrogate exchanges "
        "over that cell's own seasonal temperature range. Above it the omitted "
        "phase term dominates the sensible term the depth classes are about; "
        "below it, it does not.")
    freezing["caveat"] = (
        "the lake SET comes from surface_water.nc, which on this build is "
        "forced by the bootstrap climatology, while the temperatures are the "
        "baseline's. The all_land_area population is reported beside it so "
        "that a verdict holding on both does not rest on which cells hold "
        "lakes.")
    report["freezing_threshold"] = freezing

    # ---- LSHY-6: what the crossing loses at the consumer -----------------
    blended = (1.0 - f_lake) * capacity + f_lake * lake_dwmax_m
    beta_soil = wetness(soil_water, capacity[None, :, :])
    beta_lake = wetness(soil_water, np.full_like(capacity, lake_dwmax_m)[None, :, :])
    beta_blend = wetness(soil_water, blended[None, :, :])
    beta_tiled = ((1.0 - f_lake)[None, :, :] * beta_soil
                  + f_lake[None, :, :] * beta_lake)

    # The potential latent heat flux, CAPPED at the surface energy available.
    # hfls / beta is the exact inverse and is unusable: beta is under 0.05 on
    # nearly half this planet's land bins and the quotient reaches 1e17 W m-2.
    ceiling = np.maximum(net_radiation, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        implied = np.where(beta_soil > 0.0, latent / beta_soil, 0.0)
    potential_estimate = np.minimum(np.nan_to_num(implied), ceiling)

    def flux_gap(delta_beta: np.ndarray, mask: np.ndarray) -> dict:
        out = {}
        for label, potential in (("estimate", potential_estimate), ("upper_bound", ceiling)):
            annual = (np.nan_to_num(delta_beta * potential)
                      * bin_weight[:, None, None]).sum(axis=0)
            annual = np.where(mask, annual, 0.0)
            planet = float((annual * cell_km2).sum() / cell_km2.sum())
            out[label] = {
                "planet_mean_w_m2": planet,
                "land_mean_w_m2": float((annual * land_km2).sum() / land_km2.sum()),
                "multiples_of_the_bar": abs(planet) / BAR_W_M2,
                "worst_cell_w_m2": float(annual.flat[np.argmax(np.abs(annual))]),
            }
        return out

    # WHERE THE SIGN REVERSES. The wetness factor is flat below the saturation
    # knee at C = w / 0.4 and convex above it, so its curvature changes sign
    # AT the knee. A cell whose two tiles straddle the knee carries the
    # opposite sign to one whose tiles are both on the unsaturated branch, and
    # the gap is therefore not one-signed in the lake fraction. This counts the
    # lake area on each side rather than assuming either.
    knee = soil_water / DRHSFULL
    soil_saturated = capacity[None, :, :] <= knee
    lake_saturated = lake_dwmax_m <= knee
    straddles = (soil_saturated != lake_saturated)
    both_unsaturated = ~soil_saturated & ~lake_saturated
    lake_bins = np.broadcast_to(lake_km2[None, :, :], straddles.shape)
    weights = lake_bins * bin_weight[:, None, None]
    has_lake_bins = np.broadcast_to(has_lake[None, :, :], straddles.shape)
    report["curvature"] = {
        "finding": "the wetness factor min(1, w / (0.4 C)) is NOT convex in "
                   "the capacity. It is flat below the saturation knee at "
                   "C = w / 0.4 and convex above it, so the derivative jumps "
                   "down at the knee and the function is concave there. A "
                   "two-point capacity distribution is one-signed only when "
                   "the law it is pushed through has one sign of curvature, "
                   "and this one does not: a cell whose soil and lake tiles "
                   "straddle the knee carries the opposite sign, and the "
                   "sign can reverse within a single cell as the lake "
                   "fraction is varied. The runoff consumer, "
                   "max(0, w + F dt - C) / dt, IS convex in the capacity "
                   "throughout and is not affected.",
        "share_of_lake_area_bins_whose_tiles_straddle_the_knee": share(
            weights[has_lake_bins], straddles[has_lake_bins]),
        "share_of_lake_area_bins_with_both_tiles_unsaturated": share(
            weights[has_lake_bins], both_unsaturated[has_lake_bins]),
        "consequence": "world-cyu3's pre-registered invariant 2 holds for "
                       "runoff and fails for the wetness factor. Any tile "
                       "operator chosen for this crossing has to carry the "
                       "distribution rather than a corrected mean, because a "
                       "correction with one sign cannot track a gap with two.",
    }

    report["mosaic_cost"] = {
        "mosaic_absent_today": {
            "what": "the lake tile is not staged at all, so the cell runs on "
                    "the upland-soil capacity alone. This is the loss the "
                    "world carries today, and it is the whole tile rather "
                    "than a mixing error.",
            "all_land_cells": flux_gap(beta_tiled - beta_soil, is_land),
        },
        "mosaic_blended_if_lakes_were_staged": {
            "what": "the Jensen term that would remain: the two-point capacity "
                    "distribution averaged before the convex wetness law "
                    "rather than after it.",
            "all_land_cells": flux_gap(beta_tiled - beta_blend, is_land),
        },
        "population": "every land cell. There is no excluded population: the "
                      "flux arms stand on the capacity the climatology's own "
                      "run integrated, so the conservation test above passes "
                      "and there is nothing to report around.",
    }

    # THE DIRECTION IS A PROPERTY OF THE FIELD, so both fields are reported and
    # neither is called "the" capacity. The blend raises the capacity wherever
    # lake_dwmax_m exceeds the soil capacity and lowers it wherever it does
    # not, and the two fields do not agree about which is the common case.
    def perturbation_arm(field: np.ndarray) -> dict:
        mixed = (1.0 - f_lake) * field + f_lake * lake_dwmax_m
        with np.errstate(invalid="ignore"):
            size = np.abs(mixed - field) / np.where(field > 0, field, np.nan)
            # THE ABSENT TILE IS A MEAN-CAPACITY CHANGE and the mixing term is
            # not: the blend IS the area-weighted mean of the two-point
            # distribution, so it moves no mean at all. Only the first can be
            # taken through a chord in the mean, which is the only instrument
            # registered for the runoff consumer. world-cyu3.
            log_change = (np.log(np.where(field > 0, field, np.nan))
                          - np.log(np.where(mixed > 0, mixed, np.nan)))
        weights = land_km2[is_land]
        return {
            "land_mean_capacity_m": float(
                (field * land_km2)[is_land].sum() / land_km2[is_land].sum()),
            "absent_tile_mean_log_capacity_change": float(
                np.nansum(log_change[is_land] * weights) / weights.sum()),
            "absent_tile_mean_absolute_log_capacity_change": float(
                np.nansum(np.abs(log_change[is_land]) * weights) / weights.sum()),
            "over_lake_bearing_land_area": area_weighted_quantile(
                size[has_lake], land_km2[has_lake], [0.5, 0.75, 0.9, 0.95, 0.99]),
            "max": float(np.nanmax(size[has_lake])),
            "share_of_lake_bearing_land_area_above_100_percent": share(
                land_km2[has_lake], size[has_lake] > 1.0),
            "lake_bearing_land_area_share_where_lake_dwmax_exceeds_the_soil_capacity":
                share(land_km2[has_lake], field[has_lake] < lake_dwmax_m),
        }

    report["capacity_perturbation"] = {
        "f_lake_over_lake_bearing_cells": {
            "median": float(np.median(f_lake[has_lake])),
            "mean": float(f_lake[has_lake].mean()),
            "max": float(f_lake[has_lake].max()),
        },
        "on_the_capacity_that_ran": perturbation_arm(capacity),
        "on_the_staged_capacity": perturbation_arm(staged_capacity),
        "direction": "where lake_dwmax_m exceeds the soil capacity the blend "
                     "RAISES the cell's capacity, which raises the water "
                     "needed to reach 40 per cent of it and so lowers both the "
                     "wetness factor and saturation-excess runoff. "
                     "config/planet.yaml argues lake_dwmax_m shallower than "
                     "the capacity it is blended into, for the opposite "
                     "effect. Which way the field actually goes is measured "
                     "per field above rather than assumed, because the two "
                     "iterations do not agree. world-kvr owns the value's "
                     "derivation.",
    }

    report["provenance"] = {
        "generator": "hydrography/scripts/lake_mosaic_cost.py",
        "inputs": {
            "support": str(support_path.relative_to(PROJECT_ROOT)),
            "climatology": rel(clim_path),
            "capacity_that_ran": consumed["path"],
            "staged_capacity": staged["path"],
        },
        "input_sha256": {
            "capacity_that_ran": consumed["sha256"],
            "staged_capacity": staged["sha256"],
        },
        "capacity_pairing": consumed,
        "staged_capacity_stamp": staged,
        "reads_only": True,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=False) + "\n")
    print(f"wrote {args.out}")
    print(json.dumps(report["freezing_threshold"], indent=2)[:1200])
    print(json.dumps(report["mosaic_cost"], indent=2)[:2000])


def selftest() -> None:
    """Identities the instrument must satisfy, each with a right answer."""
    # 1. The wetness factor saturates exactly at DRHSFULL of capacity.
    assert wetness(np.array([0.4]), np.array([1.0]))[0] == 1.0
    assert abs(wetness(np.array([0.2]), np.array([1.0]))[0] - 0.5) < 1e-12

    # 2. The blend is affine, so tiling and blending AGREE at both endpoints
    #    of the lake fraction. This is world-cyu3's invariant 2 as a test.
    water = np.array([0.05])
    for fraction in (0.0, 1.0):
        soil, lake = np.array([0.09]), np.array([0.2])
        blended = (1 - fraction) * soil + fraction * lake
        tiled = (1 - fraction) * wetness(water, soil) + fraction * wetness(water, lake)
        assert abs(tiled[0] - wetness(water, blended)[0]) < 1e-12, fraction

    # 3. Where BOTH tiles sit on the unsaturated branch the law is convex in
    #    the capacity, so tiling lands at or above blending. This is the case
    #    world-cyu3's invariant 2 describes, and on this branch it holds.
    dry = np.array([0.02])          # below 0.4 * both capacities
    for fraction in (0.1, 0.25, 0.5, 0.75, 0.9):
        soil, lake = np.array([0.09]), np.array([0.2])
        blended = (1 - fraction) * soil + fraction * lake
        tiled = (1 - fraction) * wetness(dry, soil) + fraction * wetness(dry, lake)
        assert tiled[0] >= wetness(dry, blended)[0] - 1e-12, fraction

    # 4. WHERE THE TWO TILES STRADDLE THE SATURATION KNEE THE SIGN REVERSES,
    #    and the gap is not one-signed in the lake fraction at all. The wetness
    #    factor is min(1, w / (0.4 C)): flat below the knee at C = w / 0.4 and
    #    convex above it, so the derivative JUMPS DOWN at the knee and the
    #    function is CONCAVE there. A two-point distribution has one sign of
    #    curvature only when the function it is pushed through has one, and
    #    this one does not. world-cyu3's invariant 2 holds for the runoff
    #    consumer, which is max(0, ...) of a decreasing linear function of the
    #    capacity and therefore convex throughout, and fails for this one.
    straddling = np.array([0.05])   # saturates the 0.09 m tile, not the 0.2 m
    signs = []
    for fraction in (0.1, 0.25, 0.5, 0.75, 0.9):
        soil, lake = np.array([0.09]), np.array([0.2])
        blended = (1 - fraction) * soil + fraction * lake
        tiled = ((1 - fraction) * wetness(straddling, soil)
                 + fraction * wetness(straddling, lake))
        signs.append(np.sign(tiled[0] - wetness(straddling, blended)[0]))
    assert min(signs) < 0 < max(signs), signs

    # 5. Both endpoints still agree exactly on the straddling case, so the
    #    sign reversal is interior and is not a broken endpoint.
    for fraction in (0.0, 1.0):
        soil, lake = np.array([0.09]), np.array([0.2])
        blended = (1 - fraction) * soil + fraction * lake
        tiled = ((1 - fraction) * wetness(straddling, soil)
                 + fraction * wetness(straddling, lake))
        assert abs(tiled[0] - wetness(straddling, blended)[0]) < 1e-12, fraction

    # 6. A cell saturated in BOTH tiles is insensitive to capacity in both
    #    orders, so the gap is exactly zero. world-cyu3's invariant 1.
    wet = np.array([10.0])
    fraction = 0.5
    soil, lake = np.array([0.09]), np.array([0.2])
    blended = (1 - fraction) * soil + fraction * lake
    tiled = (1 - fraction) * wetness(wet, soil) + fraction * wetness(wet, lake)
    assert tiled[0] == wetness(wet, blended)[0] == 1.0

    # 7. The equivalent ice thickness inverts the audit's own figure: one
    #    metre of ice is worth a 93 K swing of the surrogate.
    swing = ICE_LATENT_HEAT_J_M2_PER_M / SURROGATE_HEAT_CAPACITY_J_M2_K
    assert abs(swing - 92.9) < 0.5, swing

    # 8. The area-weighted quantile of a constant is that constant.
    q = area_weighted_quantile(np.full(5, 7.0), np.arange(1.0, 6.0), [0.1, 0.9])
    assert all(abs(v - 7.0) < 1e-12 for v in q.values())

    print("selftest: 8 identities hold, including the sign reversal at the "
          "saturation kink")


if __name__ == "__main__":
    main()
