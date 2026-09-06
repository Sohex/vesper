#!/usr/bin/env python3
"""What a drift in the simulated cover is worth to the climate run that reads it.

    python biosphere/scripts/derive_cover_tolerance.py [--run <lpj run dir>]

Vesper is a simulated world. This derives the tolerance
`biosphere/config/equilibrium_window.yaml` puts on the end-to-end drift of the
two `fpc.out` cover quantities, by carrying that drift through the consumer that
closes a loop with it: `exoplasim/scripts/build_surface_albedo.py --mode
modelled` stages a surface albedo from the simulated foliar cover, and the next
climate run is forced with it. The tolerance is the drift whose surface
temperature is below what the climate arm has declared it will carry.

NOTHING IS DECLARED HERE. Every link is read from the artifact that owns it, so
a link that moves moves this answer, and the module's job is to say whether the
declared tolerance still agrees with the chain rather than to hold a number.

## The chain, and it is the consumer's own arithmetic

`composite_rootable` mixes cover into the surface albedo inside the rootable
area only. Below the cover ceiling -- the regime almost every gridcell is in --
its derivative with respect to a relative drift `X` applied to both covers is

    d(cell surface albedo) = X * rootable * [ tree*(a_tree - a_sub)
                                            + grass*(a_grass - a_sub) ]

so the planet-mean surface albedo moves by the area-weighted sum of that over
the gridcells BOTH models own, divided by the whole planet's area. Both, because
the ecology model simulates every gridcell with any native land and the climate
model reads a staged surface field only where its own binary mask says land:
simulated cover on a cell that mask calls ocean reaches no climate run and is
worth no kelvin. That reaches the top of the modelled atmosphere through the
MEASURED surface-to-planetary attenuation `scripts/error_budget.py` owns, and
reaches kelvin through `lib/sensitivity.py`, which is this project's one
flux-to-kelvin conversion and is not re-derived here.

## The bar, and it is read rather than chosen

`exoplasim/scripts/assess_convergence.py` declares `OFFSET_TOLERANCE_K`, the
residual offset a converged climate run is allowed to carry, and
`RESOLVING_FACTOR`, the factor by which a term must sit below a threshold before
it stops being able to flip a verdict on it. Both were fixed before they were
applied to any assessment. A drift in the simulated cover is a SECOND error term
in the same kelvin, so the share it may take of that tolerance is
`OFFSET_TOLERANCE_K / RESOLVING_FACTOR`.

Three other bars were considered and are not this one. `OFFSET_TOLERANCE_K`
itself would give a single extra term the whole of the allowance the climate arm
gives itself, which is a budget with two terms at one term's size. The standard
error of a climatology mean and the smallest difference a paired A/B off a
common donor resolves are NOISE FLOORS -- what the arm can see -- and a bias
below the noise floor is a far stricter demand than a bias that changes no
verdict; the paired one is also about a difference in which a drift common to
both arms cancels exactly, so it is a bar on a quantity this drift does not
enter.

## What the answer can and cannot resolve

The chain's own bracket spans a factor of about two, so it cannot distinguish
tolerances closer together than that, and a verdict finer than that would be
believing a number smaller than its own instrument's scatter. `AGREEMENT_FACTOR`
is the factor within which the declared tolerance is retained and it was
registered before the chain was taken with a measured attenuation.

The cover magnitudes come from an LPJ run, and no run in this tree is accepted
yet, so `cover_robustness` reports how far those magnitudes would have to move
to change the verdict. That is what makes an answer taken on a refused run worth
having: the verdict, not the covers, is what has to survive.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "lib"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT / "exoplasim" / "scripts"))

import numpy as np  # noqa: E402
import yaml  # noqa: E402

import builds  # noqa: E402
import gridding  # noqa: E402
import rungs  # noqa: E402
import sensitivity  # noqa: E402

CONFIG = PROJECT_ROOT / "config" / "planet.yaml"
CONTRACT = PROJECT_ROOT / "biosphere" / "config" / "equilibrium_window.yaml"

# The factor within which the declared tolerance is retained rather than moved.
# REGISTERED before the chain was taken against a measured attenuation, and it
# is not a preference: a chain of this many bracketed links cannot resolve
# better than this, so a verdict inside it would be reading the instrument's own
# scatter as a result. `resolving_factor_of_the_chain` on the report is what the
# bracket actually spans, and it is reported so this number stays checkable.
AGREEMENT_FACTOR = 2.0

# The two quantities the chain is about, by their contract ids.
COVER_QUANTITIES = ("fpc.out tree cover", "fpc.out grass cover")

GRASS_COLUMNS = ("C3G", "C4G")
NOT_A_TYPE = ("Lon", "Lat", "Year", "Total")


def climate_bar_k() -> tuple[float, dict]:
    """The share of the climate arm's own offset tolerance one term may take."""
    import assess_convergence

    tolerance = float(assess_convergence._OFFSET_TOLERANCE_K)
    factor = float(assess_convergence._RESOLVING_FACTOR)
    return tolerance / factor, {
        "offset_tolerance_k": tolerance,
        "resolving_factor": factor,
        "source": "exoplasim/scripts/assess_convergence.py",
        "meaning": ("the residual offset a converged climate run may carry, "
                    "divided by the factor at which a term stops being able to "
                    "flip a verdict on it"),
    }


def attenuation() -> tuple[float, tuple[float, float], dict]:
    """The MEASURED surface-to-planetary attenuation and its measured span."""
    import error_budget

    value = float(error_budget.DEFAULT_ATTENUATION)
    lo, hi = (float(x) for x in error_budget.ATTENUATION_SPREAD)
    return value, (lo, hi), {
        "source": "scripts/error_budget.py",
        "finding": "notes/audits/albedo-attenuation.md",
    }


# How far this chain's model-land area fraction may sit from the one the
# attenuation was measured against before the two are not the same planet. The
# attenuation rows were back-solved on that fraction, so a chain taken over a
# different land mask would be multiplying a kelvin by an attenuation derived
# for somewhere else. Four decimals is finer than any bracket here and coarser
# than the float the two arrive as.
LAND_FRACTION_TOLERANCE = 1.0e-4


def check_land_fraction(measured: float) -> None:
    """The one link this chain shares with a number computed somewhere else.

    `scripts/error_budget.py:ATTENUATION_LAND_FRACTION` is the model-land area
    fraction the measured attenuation's own arms were taken over, computed there
    from the staged fields rather than from this mask. Recomputing it here off
    the land mask and comparing is a check with a right answer: the two are the
    same quantity, and a disagreement says the chain and the attenuation it
    multiplies belong to different builds.
    """
    import error_budget

    expected = float(error_budget.ATTENUATION_LAND_FRACTION)
    if abs(measured - expected) > LAND_FRACTION_TOLERANCE:
        raise SystemExit(
            f"this chain's model-land area fraction {measured:.6f} disagrees "
            f"with the {expected:.6f} the measured attenuation was taken over "
            "(scripts/error_budget.py:ATTENUATION_LAND_FRACTION). The chain and "
            "the attenuation it multiplies are not on the same land mask.")


def cover_means(run_dir: Path, nlat: int, nlon: int):
    """Area-weighted planet-mean tree and grass cover, per modelled gridcell.

    Returned as PLANET means rather than land means, so the chain needs no
    separate land-fraction link: a gridcell the run does not simulate
    contributes nothing to a cover change because it carries no simulated cover.

    Cells are matched to the model grid through `gridding.model_label_cells`,
    which refuses a label that is not on the model's own axis rather than
    snapping it to the nearest column. The longitudes an LPJ table writes are
    signed, so they are taken modulo 360 first; that is an exact identity on the
    label and not a nearest match.
    """
    import pandas as pd

    path = Path(run_dir) / "fpc.out"
    if not path.is_file():
        raise SystemExit(f"{path} does not exist")
    weights = gridding.gaussian_grid(nlat, nlon).cell_area_fraction()
    tree_sum = np.zeros((nlat, nlon))
    grass_sum = np.zeros((nlat, nlon))
    counts = np.zeros((nlat, nlon))
    columns: list[str] | None = None
    for chunk in pd.read_csv(path, sep=r"\s+", chunksize=400_000):
        if columns is None:
            columns = list(chunk.columns)
            tree_cols = [c for c in columns
                         if c not in NOT_A_TYPE and c not in GRASS_COLUMNS]
            grass_cols = [c for c in GRASS_COLUMNS if c in columns]
            if not tree_cols or not grass_cols:
                raise SystemExit(f"{path} has no tree or no grass columns")
        rows, cols = gridding.model_label_cells(
            chunk["Lat"].to_numpy(), chunk["Lon"].to_numpy() % 360.0,
            nlat, nlon, what=f"a row of {path.name}")
        flat = rows * nlon + cols
        np.add.at(tree_sum.reshape(-1), flat,
                  chunk[tree_cols].to_numpy().sum(axis=1))
        np.add.at(grass_sum.reshape(-1), flat,
                  chunk[grass_cols].to_numpy().sum(axis=1))
        np.add.at(counts.reshape(-1), flat, 1.0)
    simulated = counts > 0
    if not simulated.any():
        raise SystemExit(f"{path} matched no model gridcell")
    tree = np.where(simulated, tree_sum / np.maximum(counts, 1.0), 0.0)
    grass = np.where(simulated, grass_sum / np.maximum(counts, 1.0), 0.0)
    return tree, grass, simulated, weights


def rootable_field(config: dict, rung: str, nlat: int, nlon: int,
                   simulated: np.ndarray) -> np.ndarray:
    """BIO-11's rootable fraction, and it is REQUIRED wherever the run simulated.

    The artifact carries no rootable fraction where the model grid has no native
    land, so the field arrives masked. A masked cell the run never simulated
    contributes nothing and is filled with zero; a masked cell the run DID
    simulate is a disagreement between the rootable partition and the driver
    about where this world's land is, and it refuses rather than being filled --
    filling it would quietly leave the cover of a real gridcell out of the
    planet mean this chain is taken over.
    """
    from netCDF4 import Dataset

    target = (builds.component_data("biosphere", config)
              / f"rootable_fraction_{rung}.nc")
    if not target.is_file():
        raise SystemExit(
            f"{target} does not exist. BIO-11's rootable fraction is a link in "
            "this chain and is not substituted for.")
    with Dataset(target) as ds:
        field = np.asarray(ds["f_rootable"][:], dtype=float)
    if field.shape != (nlat, nlon):
        raise SystemExit(f"{target} is not on the {nlat}x{nlon} model grid")
    absent = ~np.isfinite(field)
    if np.any(absent & simulated):
        raise SystemExit(
            f"{target} carries no rootable fraction at "
            f"{int((absent & simulated).sum())} gridcells the run simulated")
    return np.where(absent, 0.0, field)


def model_land(rung: str, nlat: int, nlon: int) -> np.ndarray:
    """The cells the climate model itself calls land, from the staged mask.

    THIS IS THE SET THE CHAIN IS TAKEN OVER AND IT IS SMALLER THAN THE SET THE
    ECOLOGY MODEL RUNS ON. LPJ simulates every gridcell with any native land;
    ExoPlaSim reads a staged surface field only where its own binary mask says
    land, so simulated cover on a cell the mask calls ocean reaches no climate
    run and is worth no kelvin. Weighting the chain over the simulated cells
    instead inflates it by better than a factor of two, which is larger than
    every bracket in the chain put together.
    """
    sys.path.insert(0, str(PROJECT_ROOT / "exoplasim" / "scripts"))
    import sra

    target = (PROJECT_ROOT / "exoplasim" / "inputs" / rung.lower()
              / f"orogen_{rung}_surf_0172.sra")
    if not target.is_file():
        raise SystemExit(
            f"{target} does not exist. The model's own land mask decides which "
            "simulated cover reaches a climate run and is not reconstructed.")
    mask = sra.read_sra(target, nlat, nlon)
    return np.asarray(mask) > 0.5


def kelvin_per_unit_drift(tree, grass, rootable, simulated, weights,
                          a_tree: float, a_grass: float, a_sub: float,
                          att: float, alpha: float, slope: float) -> float:
    """Kelvin per unit RELATIVE drift applied to both simulated covers."""
    per_cell = rootable * (tree * (a_tree - a_sub) + grass * (a_grass - a_sub))
    d_surface = float((per_cell * weights * simulated).sum() / weights.sum())
    d_toa = d_surface * att
    return abs(sensitivity.planetary_albedo_to_kelvin(d_toa, alpha,
                                                      slope=slope))


def derive(run_dir: Path, config: dict | None = None) -> dict:
    config = config or yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    contract = yaml.safe_load(CONTRACT.read_text(encoding="utf-8"))
    rung, nlat, nlon = rungs.model_grid(config)

    model = config["model"]
    a_tree = float(model["tree_albedo"])
    a_grass = float(model["grass_albedo"])
    tree_bracket = tuple(float(x) for x in model["tree_albedo_bracket"])
    grass_bracket = tuple(float(x) for x in model["grass_albedo_bracket"])

    report_path = (PROJECT_ROOT / "exoplasim" / "inputs"
                   / rung.lower() / "albedo_report.json")
    if not report_path.is_file():
        raise SystemExit(f"{report_path} does not exist; the substrate albedo "
                         "this chain differences against is not substituted for")
    a_sub = float(json.loads(report_path.read_text(encoding="utf-8"))
                  ["land_mean_bare_rock"])

    tree, grass, ran, weights = cover_means(run_dir, nlat, nlon)
    land = model_land(rung, nlat, nlon)
    # The cells whose simulated cover reaches a climate run: both models' own.
    simulated = ran & land
    rootable = rootable_field(config, rung, nlat, nlon, simulated)
    att, att_span, att_source = attenuation()
    alpha, alpha_source = sensitivity.planetary_albedo(config)
    slope = float(sensitivity.SLOPE_K_PER_FLUX_RATIO)
    slope_span = tuple(float(x) for x in sensitivity.SLOPE_SPREAD_K_PER_FLUX_RATIO)
    bar, bar_source = climate_bar_k()

    land_fraction = float((weights * land).sum() / weights.sum())
    check_land_fraction(land_fraction)

    def k(att_, at, ag, sl):
        return kelvin_per_unit_drift(tree, grass, rootable, simulated, weights,
                                     at, ag, a_sub, att_, alpha, sl)

    central = k(att, a_tree, a_grass, slope)
    # Two identities the chain has to satisfy, checked rather than assumed. Both
    # catch a term that has stopped being proportional to what it is written as.
    doubled = kelvin_per_unit_drift(2.0 * tree, 2.0 * grass, rootable,
                                    simulated, weights, a_tree, a_grass, a_sub,
                                    att, alpha, slope)
    if not math.isclose(doubled, 2.0 * central, rel_tol=1e-9):
        raise SystemExit("the chain is not linear in the simulated cover; its "
                         "derivative is not what the consumer's own formula has")
    if not math.isclose(k(2.0 * att, a_tree, a_grass, slope), 2.0 * central,
                        rel_tol=1e-9):
        raise SystemExit("the chain is not linear in the attenuation; a "
                         "measured span cannot be carried through it")
    corners = [k(x, at, ag, sl)
               for x in att_span for at in tree_bracket
               for ag in grass_bracket for sl in slope_span]
    k_lo, k_hi = min(corners), max(corners)
    if k_lo <= 0.0:
        raise SystemExit("the chain returns no temperature response; the "
                         "simulated cover or the albedo contrast is zero")

    declared = sorted({float(q["relative_end_to_end_limit"])
                       for q in contract["assessed"]["quantities"]
                       if q["id"] in COVER_QUANTITIES})
    if len(declared) != 1:
        raise SystemExit(
            "the two cover quantities carry different tolerances "
            f"{declared}; this chain derives one number for both")
    declared_limit = declared[0]

    derived_limit = bar / central
    derived_span = (bar / k_hi, bar / k_lo)
    ratio = declared_limit / derived_limit
    agrees = 1.0 / AGREEMENT_FACTOR <= ratio <= AGREEMENT_FACTOR

    # How far the simulated covers would have to move to change the verdict.
    # The response is linear in them, so the bound is exact in the ratio.
    cover_up = AGREEMENT_FACTOR / ratio
    cover_down = 1.0 / (AGREEMENT_FACTOR * ratio)

    return {
        "generator": "biosphere/scripts/derive_cover_tolerance.py",
        "run": str(run_dir),
        "contract_version": contract["contract_version"],
        "quantities": list(COVER_QUANTITIES),
        "links": {
            "tree_albedo": [a_tree, list(tree_bracket), "config/planet.yaml"],
            "grass_albedo": [a_grass, list(grass_bracket), "config/planet.yaml"],
            "substrate_albedo": [a_sub, None,
                                 f"exoplasim/inputs/{rung.lower()}"
                                 "/albedo_report.json land_mean_bare_rock"],
            "rootable_fraction_land_mean": [
                float((rootable * weights * simulated).sum()
                      / (weights * simulated).sum()), None,
                "biosphere/data/<build>/rootable_fraction_<rung>.nc"],
            "planet_mean_tree_cover": [
                float((tree * weights * simulated).sum() / weights.sum()), None,
                f"{run_dir}/fpc.out"],
            "planet_mean_grass_cover": [
                float((grass * weights * simulated).sum() / weights.sum()), None,
                f"{run_dir}/fpc.out"],
            "attenuation": [att, list(att_span), att_source["source"]],
            "flux_to_kelvin": [slope, list(slope_span), "lib/sensitivity.py"],
            "planetary_albedo": [alpha, None, alpha_source["source"]],
        },
        "cells": {
            "simulated_by_the_ecology_model": int(ran.sum()),
            "called_land_by_the_climate_model": int(land.sum()),
            "both, and the chain is taken over these": int(simulated.sum()),
            "climate_land_the_run_did_not_simulate": int((land & ~ran).sum()),
            "model_land_area_fraction": land_fraction,
        },
        "bar_k": bar,
        "bar": bar_source,
        "kelvin_per_unit_relative_drift": central,
        "kelvin_per_unit_relative_drift_span": [k_lo, k_hi],
        "declared_limit": declared_limit,
        "declared_limit_kelvin": declared_limit * central,
        "declared_limit_kelvin_span": [declared_limit * k_lo,
                                       declared_limit * k_hi],
        "derived_limit": derived_limit,
        "derived_limit_span": list(derived_span),
        "resolving_factor_of_the_chain": derived_span[1] / derived_span[0],
        "agreement_factor": AGREEMENT_FACTOR,
        "declared_over_derived": ratio,
        "verdict": "retained" if agrees else "outside the registered factor",
        "cover_robustness": {
            "covers_larger_by": cover_up,
            "covers_smaller_by": cover_down,
            "meaning": ("the factor the simulated covers would have to move by "
                        "for the declared tolerance to leave the registered "
                        "agreement factor; the response is linear in them"),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run", type=Path, required=True,
                        help="an LPJ run directory holding fpc.out")
    parser.add_argument("--json", action="store_true",
                        help="print the whole report rather than the summary")
    args = parser.parse_args()

    result = derive(args.run)
    if args.json:
        print(json.dumps(result, indent=2))
        return 0 if result["verdict"] == "retained" else 1

    k = result["kelvin_per_unit_relative_drift"]
    lo, hi = result["kelvin_per_unit_relative_drift_span"]
    print(f"chain: {k:.4f} K per unit relative cover drift "
          f"[{lo:.4f}, {hi:.4f}]")
    print(f"the declared {result['declared_limit']} is worth "
          f"{result['declared_limit_kelvin']:.4f} K "
          f"[{result['declared_limit_kelvin_span'][0]:.4f}, "
          f"{result['declared_limit_kelvin_span'][1]:.4f}]")
    print(f"the bar is {result['bar_k']:.4f} K "
          f"({result['bar']['offset_tolerance_k']} / "
          f"{result['bar']['resolving_factor']}, "
          f"{result['bar']['source']})")
    print(f"derived limit {result['derived_limit']:.4f} "
          f"[{result['derived_limit_span'][0]:.4f}, "
          f"{result['derived_limit_span'][1]:.4f}], a factor of "
          f"{result['resolving_factor_of_the_chain']:.2f}")
    print(f"declared / derived = {result['declared_over_derived']:.3f} against "
          f"a registered factor of {result['agreement_factor']}: "
          f"{result['verdict'].upper()}")
    rob = result["cover_robustness"]
    print(f"the simulated covers would have to be "
          f"{rob['covers_larger_by']:.2f}x larger or "
          f"{rob['covers_smaller_by']:.2f}x smaller to change that")
    return 0 if result["verdict"] == "retained" else 1


if __name__ == "__main__":
    raise SystemExit(main())
