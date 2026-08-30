#!/usr/bin/env python3
"""Is a per-cell land emissivity worth carrying, or is a scalar enough?

    python analysis/emissivity_contrast.py --land-mean   # the scalar
    python analysis/emissivity_contrast.py                # is a field worth it

Worldbuilding. Vesper is an invented planet; every number below is a property
of the simulated land surface or of the model that radiates from it.

`analysis/rock_emissivity.py` derives a broadband thermal emissivity per Orogen
rock class. This decides whether that variation is large enough to be worth
giving the model as a per-cell surface field, rather than one scalar in
`radmod_nl`. It is a test with a right answer rather than a comparison: the
criterion is fixed below, in the units the model's own surface energy balance
reports, and the script says pass or fail against it.

## The quantity, and why it is exactly computable

`radmod.f90:lwr` makes the surface net longwave

    net = -eps * (sigma*Ts^4 - LWdown)

which is exactly linear in the emissivity and in nothing else the surface
carries. So a contrast of `d` in emissivity is worth `d` times the model's own
surface longwave loss, output code 177, per cell -- no scheme, no feedback, no
run. That is what makes the size of the effect knowable in advance, which is
what CLAUDE.md's rule about checking the instrument against the size of the
effect asks for.

Two separate quantities come out of that, and conflating them is how a field
gets built for the wrong reason:

**The mean.** What one scalar at the lithology-weighted land mean is worth
against whatever scalar is declared. This is a property of the NUMBER, and a
field is not needed to fix it. `--land-mean` writes that scalar and its terrain
hash to `analysis/land_emissivity.json` and reads NO climatology, because the
mean is an area weighting over the mesh's own land and nothing else.

**The contrast.** What the field is worth OVER a scalar already set to the
field's own land mean. This is the only thing a per-cell field buys, and it is
what the criterion below is applied to.

## The criterion, fixed before any spectrum was read

The contrast must exceed 1.4 W/m2 in surface net longwave. That is the top of
the model's own dry adiabatic energy sink, which the filter hides and the
surface fluxes pay for: a surface energy term smaller than it is inside the
model's own non-conservation and cannot be claimed as a resolved effect however
tidy its derivation looks.

It is applied to the area-weighted standard deviation over land, not to the
range, because a range is set by two cells and a field is bought for all of
them. The share of land above the threshold is reported beside it, because a
term that fails on the mean and passes on a small area is a different finding
from one that fails everywhere -- and this world's closed-basin floors are
exactly the small area that has earned a field before.

## The bound that can overturn the answer

The spectra stop near 14 micrometres and a land surface at these temperatures
radiates most of its energy beyond that. `rock_emissivity.py` carries the
in-band value to the whole thermal spectrum. The opposing bound is a blackbody
outside the band, which multiplies every between-class contrast by the in-band
Planck fraction and is therefore the bound that COMPRESSES the signal. A
verdict that holds at the optimistic bound and fails at the compressing one is
not a verdict, so both are computed and the answer is the weaker.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml
from netCDF4 import Dataset

import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))

from builds import grid_export, mesh_export  # noqa: E402
from climatology import annual_mean  # noqa: E402
from gridding import cell_fraction, cell_mean, region_cells  # noqa: E402
from orogen import Export, LAND  # noqa: E402
from paths import best_available_climatology, rel  # noqa: E402

CONFIG = ROOT / "config" / "planet.yaml"
DERIVATION = ROOT / "analysis" / "rock_emissivity.json"
REPORT = ROOT / "analysis" / "emissivity_contrast.json"
LAND_MEAN_REPORT = ROOT / "analysis" / "land_emissivity.json"

# THE CRITERION. Fixed before any spectrum was read. See the module docstring:
# the model's own dry adiabatic energy sink is 0.6 to 1.4 W/m2 and the surface
# fluxes are what pays for it, so this is the top of that range.
THRESHOLD_W_M2 = 1.4


def class_emissivity(mesh: Export, values: dict) -> np.ndarray:
    """The per-class emissivity mapped onto the mesh's cells, NaN off the table."""
    rock = mesh.field("substrate_class").astype(int)
    codes = [r["code"] for r in mesh.manifest["lithology"]["rockClasses"]]
    out = np.full(rock.shape, np.nan)
    for index, code in enumerate(codes):
        if code in values:
            out[rock == index] = values[code]
    # RULE 1: land is surface_class, never land_mask.
    is_land = mesh.surface_class == LAND
    if np.isnan(out[is_land]).any():
        missing = sorted({codes[i] for i in np.unique(rock[is_land])} - set(values))
        raise SystemExit(f"no emissivity for rock class(es) {missing}; add them "
                         "to SELECTORS in analysis/rock_emissivity.py")
    return out


def write_land_mean(args, config: dict, derivation: dict) -> None:
    """The scalar `model.land_longwave_emissivity` IS, and nothing else.

    NO CLIMATOLOGY IS INVOLVED and that is worth saying, because the contrast
    computation below needs one and the two answers are easy to confuse. The
    land mean is the area weighting of the per-class emissivity over the mesh's
    own land cells: a property of the mesh and of `rock_emissivity.json` alone,
    so it is re-derivable on any tree that has the build, before any run exists.
    What needs a climatology is the question of whether the FIELD is worth
    carrying over that scalar, because that is priced in the model's own surface
    longwave loss.

    The report carries the terrain hash because this number is DEFINED as a
    weighting over one build's land, so every terrain change moves it. A check
    that compared only the number would pass an artifact computed on a build the
    tree no longer points at, which is how the configured scalar came to carry
    another terrain's mean.
    """
    mesh = Export(args.mesh or mesh_export(config))
    values = {code: float(entry["emissivity"])
              for code, entry in derivation["classes"].items()}
    is_land = mesh.surface_class == LAND
    area = mesh.cell_area.astype(np.float64)[is_land]
    eps = class_emissivity(mesh, values)[is_land]

    def weighted(table: dict) -> float:
        per_cell = class_emissivity(
            mesh, {code: float(v) for code, v in table.items()})[is_land]
        return float((per_cell * area).sum() / area.sum())

    # The preparation arms, which are the bracket to run rather than a guessed
    # one: every class at its solid samples, and every class at its particulate
    # ones. Read from the same derivation so the bracket moves with it.
    arms, arm_gaps = {}, {}
    for arm in ("solid", "particulate"):
        table, missing = {}, []
        for code, entry in derivation["classes"].items():
            side = entry.get(arm) or {}
            if isinstance(side.get("mean"), (int, float)):
                table[code] = side["mean"]
            else:
                # A class with no samples of that preparation keeps its adopted
                # value: the arm is a bound on the PREPARATION and a class the
                # library has never measured that way cannot be moved by it.
                table[code] = entry["emissivity"]
                missing.append(code)
        arms[arm] = round(weighted(table), 4)
        if missing:
            arm_gaps[arm] = sorted(missing)

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": "analysis/emissivity_contrast.py --land-mean",
        "quantity": "model.land_longwave_emissivity: the area-weighted mean of "
                    "the per-rock-class thermal emissivity over this build's "
                    "land, on the mesh. No climatology enters it.",
        "terrain_hash": mesh.terrain_hash,
        "source_build": config["source_build"],
        "derivation": rel(args.derivation),
        "land_mean": round(float((eps * area).sum() / area.sum()), 6),
        "preparation_arms": arms,
        "preparation_arm_classes_without_samples": arm_gaps,
        "land_cells": int(is_land.sum()),
    }
    args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"land-mean emissivity {report['land_mean']:.6f} on "
          f"{config['source_build']} ({mesh.terrain_hash[:8]}), "
          f"{report['land_cells']} land cells")
    if arms:
        print(f"  preparation arms {arms}")
    print(f"wrote {args.json}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--derivation", type=Path, default=DERIVATION)
    ap.add_argument("--climatology", type=Path, default=None,
                    help="a climatology carrying rls (surface net longwave, "
                         "code 177) and lsm, on the grid the config names. "
                         "Defaults to the best available, which is the "
                         "baseline when config names one and the bootstrap "
                         "when it does not. Not read at all by --land-mean")
    ap.add_argument("--land-mean", action="store_true",
                    help="write only the land-mean scalar, which needs no "
                         "climatology, and stop")
    ap.add_argument("--mesh", type=Path, default=None)
    ap.add_argument("--grid", type=Path, default=None)
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()
    args.json = args.json or (LAND_MEAN_REPORT if args.land_mean else REPORT)

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    if args.land_mean:
        write_land_mean(args, config,
                        json.loads(args.derivation.read_text(encoding="utf-8")))
        return
    # THE BEST AVAILABLE, not the first available. What a per-cell field is
    # worth is priced in the model's own surface longwave loss, so it is a
    # property of one run's climate state and there is no Earth fallback for it;
    # the resolver returns the baseline when the configuration names one and the
    # bootstrap when it does not, and says which, so the report records the
    # stage the verdict rests on rather than leaving it to the caller's memory.
    args.climatology, clim_stage = best_available_climatology(args.climatology)
    if not args.climatology.is_file():
        raise SystemExit(
            f"{args.climatology} is not on disk. The contrast needs a "
            "climatology carrying rls and lsm; --land-mean writes the scalar "
            "and needs none.")
    derivation = json.loads(args.derivation.read_text(encoding="utf-8"))
    values = {code: float(entry["emissivity"])
              for code, entry in derivation["classes"].items()}

    mesh = Export(args.mesh or mesh_export(config))
    grid_dir = args.grid or grid_export(config)

    # RULE 1: land is surface_class, never land_mask.
    is_land = mesh.surface_class == LAND
    area = mesh.cell_area.astype(np.float64)
    endorheic = mesh.field("is_endorheic").astype(bool)
    eps_region = class_emissivity(mesh, values)

    # RULE 3: the mesh reaches the grid by INDEX. No longitude is compared.
    cells, nlat, nlon = region_cells(mesh, grid_dir)
    n = nlat * nlon
    mean, land_area, covered = cell_mean(cells, n, area, eps_region, is_land)
    endo, _ = cell_fraction(cells, n, area, endorheic, is_land)
    field = mean.reshape(nlat, nlon)
    endo = endo.reshape(nlat, nlon)
    land_cells = covered.reshape(nlat, nlon)

    with Dataset(args.climatology) as ds:
        rls = np.asarray(ds["rls"][:])
        lsm = np.asarray(ds["lsm"][:])
        lat = np.asarray(ds["lat"][:])
        centres = (np.asarray(ds["time"][:], dtype=float)
                   if "time" in ds.variables else None)
        clim_label = args.climatology.name
    if rls.ndim == 3 and centres is None:
        raise SystemExit(
            f"{args.climatology} carries a time axis on `rls` and no `time` "
            "variable to weight it by. The bins are not equal length, so there "
            "is no annual mean to take here; nothing reconstructs one.")
    # The time axis is reduced by the raw records each bin holds, not by one
    # over the bin count: pyburn's bins are not equal length. `rls` enters
    # linearly here, so this is a weighting correction and not a Jensen one --
    # but `lsm > 0.5` below is a threshold, and a threshold on a mis-weighted
    # mean is not a threshold on the mean.
    if rls.ndim == 3:
        rls = annual_mean(rls, centres)
    if lsm.ndim == 3:
        lsm = annual_mean(lsm, centres)
    if rls.shape != (nlat, nlon):
        raise SystemExit(
            f"{args.climatology} is {rls.shape}, the export grid is "
            f"{(nlat, nlon)}. The two must be one grid; nothing here "
            "reconstructs a coordinate.")

    # The model's own surface longwave loss, which is what an emissivity
    # contrast multiplies. rls is signed downward, so the loss is its negative.
    loss = -rls
    # Both masks are the model's: the climatology's lsm is what the run
    # integrated, and the export's land is what the emissivity is defined on.
    land = land_cells & (lsm > 0.5)
    weight = np.cos(np.deg2rad(lat))[:, None] * np.ones_like(rls) * land

    def wmean(x):
        return float((x * weight).sum() / weight.sum())

    eps_mean = wmean(field)
    loss_mean = wmean(loss)
    declared = float(config["model"].get("land_longwave_emissivity", 1.0))

    mean_term = wmean((declared - field) * loss)
    contrast = (eps_mean - field) * loss
    sd = float(np.sqrt(wmean((contrast - wmean(contrast)) ** 2)))

    def share_above(thr):
        return float(weight[np.abs(contrast) > thr].sum() / weight.sum())

    in_band = float(derivation["planck_fraction_in_band"])
    verdict = {
        "carried_to_full_spectrum": sd >= THRESHOLD_W_M2,
        "blackbody_beyond_band": sd * in_band >= THRESHOLD_W_M2,
    }
    verdict["field_is_worth_carrying"] = all(verdict.values())

    endo_heavy = land & (endo >= 0.5)
    w_endo = np.cos(np.deg2rad(lat))[:, None] * np.ones_like(rls) * endo_heavy

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "question": "does a per-cell land emissivity buy more, over a scalar "
                    "set at its own land mean, than the model's surface energy "
                    "balance can distinguish",
        "criterion_w_m2": THRESHOLD_W_M2,
        "criterion": "the area-weighted standard deviation of the surface net "
                     "longwave contrast over land, at BOTH out-of-band bounds",
        "terrain_hash": mesh.terrain_hash,
        "grid": [int(nlat), int(nlon)],
        "climatology": clim_label,
        "climatology_stage": clim_stage,
        "surface_longwave_loss_land_mean_w_m2": round(loss_mean, 2),
        "emissivity": {
            "lithology_land_mean": round(eps_mean, 4),
            "declared_scalar": declared,
            "grid_min": round(float(field[land].min()), 4),
            "grid_max": round(float(field[land].max()), 4),
        },
        "mean_term_w_m2": round(mean_term, 2),
        "contrast_w_m2": {
            "sd": round(sd, 3),
            "sd_blackbody_beyond_band": round(sd * in_band, 3),
            "max_abs": round(float(np.abs(contrast[land]).max()), 3),
            "land_area_share_above_criterion": round(share_above(THRESHOLD_W_M2), 4),
            "endorheic_heavy_mean_abs":
                round(float((np.abs(contrast) * w_endo).sum() / w_endo.sum()), 3)
                if w_endo.sum() else None,
            "all_land_mean_abs": round(wmean(np.abs(contrast)), 3),
        },
        "verdict": verdict,
    }
    args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"surface longwave loss, land mean {loss_mean:.1f} W/m2 "
          f"({clim_label})")
    print(f"scalar term: declared {declared} against a lithology land mean of "
          f"{eps_mean:.4f} is worth {mean_term:.2f} W/m2")
    print(f"field contrast over that scalar: sd {sd:.2f} W/m2, "
          f"{sd * in_band:.2f} at the compressing bound, max "
          f"{np.abs(contrast[land]).max():.2f}")
    print(f"  {share_above(THRESHOLD_W_M2) * 100:.1f}% of land above the "
          f"{THRESHOLD_W_M2} W/m2 criterion")
    print(f"verdict: a per-cell field "
          f"{'IS' if verdict['field_is_worth_carrying'] else 'is NOT'} "
          f"worth carrying at this grid")
    print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
