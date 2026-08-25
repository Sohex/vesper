#!/usr/bin/env python3
"""Is a per-cell land emissivity worth carrying, or is a scalar enough?

    python analysis/emissivity_contrast.py --climatology <file>

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
field is not needed to fix it.

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
from gridding import cell_fraction, cell_mean, region_cells  # noqa: E402
from orogen import Export, LAND  # noqa: E402

CONFIG = ROOT / "config" / "planet.yaml"
DERIVATION = ROOT / "analysis" / "rock_emissivity.json"
REPORT = ROOT / "analysis" / "emissivity_contrast.json"

# THE CRITERION. Fixed before any spectrum was read. See the module docstring:
# the model's own dry adiabatic energy sink is 0.6 to 1.4 W/m2 and the surface
# fluxes are what pays for it, so this is the top of that range.
THRESHOLD_W_M2 = 1.4


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--derivation", type=Path, default=DERIVATION)
    ap.add_argument("--climatology", type=Path, required=True,
                    help="a climatology carrying rls (surface net longwave, "
                         "code 177) and lsm, on the grid the config names. "
                         "Required and never defaulted: the loss scale is a "
                         "property of one run on one terrain")
    ap.add_argument("--mesh", type=Path, default=None)
    ap.add_argument("--grid", type=Path, default=None)
    ap.add_argument("--json", type=Path, default=REPORT)
    args = ap.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    derivation = json.loads(args.derivation.read_text(encoding="utf-8"))
    values = {code: float(entry["emissivity"])
              for code, entry in derivation["classes"].items()}

    mesh = Export(args.mesh or mesh_export(config))
    grid_dir = args.grid or grid_export(config)

    # RULE 1: land is surface_class, never land_mask.
    is_land = mesh.surface_class == LAND
    area = mesh.cell_area.astype(np.float64)
    rock = mesh.field("substrate_class").astype(int)
    endorheic = mesh.field("is_endorheic").astype(bool)
    codes = [r["code"] for r in mesh.manifest["lithology"]["rockClasses"]]

    eps_region = np.full(rock.shape, np.nan)
    for index, code in enumerate(codes):
        if code in values:
            eps_region[rock == index] = values[code]
    if np.isnan(eps_region[is_land]).any():
        missing = sorted({codes[i] for i in np.unique(rock[is_land])} - set(values))
        raise SystemExit(f"no emissivity for rock class(es) {missing}; add them "
                         "to SELECTORS in analysis/rock_emissivity.py")

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
        clim_label = args.climatology.name
    if rls.ndim == 3:
        rls = rls.mean(axis=0)
    if lsm.ndim == 3:
        lsm = lsm.mean(axis=0)
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
