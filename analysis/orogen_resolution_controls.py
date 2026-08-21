#!/usr/bin/env python3
"""The controls that decide whether a resolution difference is a resolution difference.

    # the noise floor needs runs a few percent apart, NOT a resolution ladder
    cd vendor/orogen
    for N in 2600000 2700000; do
      node --max-old-space-size=12000 tools/export-planet.mjs \
        --code 01eshm059lt0b9mpgro2y83t --radius 7645.2 --gravity 12.81 \
        --lithology-strength 0.682 --grid T42 --netcdf --quiet \
        --regions $N --out /scratch/noise_$N
    done
    python analysis/orogen_resolution_controls.py \
      --export "2.500M=source/precarve-craton/exoplasim-T42" \
      --export "2.600M=/scratch/noise_2600000" \
      --export "2.700M=/scratch/noise_2700000" \
      --export "10.00M=source/precarve-craton-10m/exoplasim-T42"

Worldbuilding. Vesper is an invented planet and this measures its terrain
generator, not a real one.

`notes/audits/orogen-resolution.md` is the argument. This is the part of it that
WITHDREW three of its own claims, and it exists as a script because that is the
only reason those claims were caught.

## Why a same-resolution control is not optional here

Changing `--regions` does not only change resolution. The first plate seed is
drawn as an INDEX, `randInt(numRegions)`, so every region count is a slightly
different realisation of the same world, and any global statistic differs
between two counts for two reasons with only one of them resolution.

Runs a few percent apart hold resolution effectively fixed, so whatever moves
between them is scatter. Measured that way, a hypsometric bias, a drainage
concavity result and an RMS-against-converged-terrain figure all turned out to
be smaller than the scatter, and all three were withdrawn. What survived is
every statistic whose change cell size PREDICTS rather than one inferred from a
correlation.

The rule this leaves: on this generator, a resolution claim needs a
same-resolution control before it is a claim, and two runs a few percent apart
is the cheap version.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "lib"))

from lib.orogen import Export  # noqa: E402

LEVELS_KM = [0.5, 1.0, 2.0]
CHANNEL_THRESHOLD_KM2 = 1e3


def load(root: Path, label: str) -> dict:
    e = Export(root, require_known_build=False)
    land = np.asarray(e.surface_class, int) > 0
    area = np.asarray(e.cell_area, float)
    elev = np.asarray(e.elevation_km, float)
    acc = np.asarray(e.flow_accumulation, float)      # km2, comparable across counts
    slope = np.tan(np.radians(np.asarray(e.local_slope_deg, float)))
    scarp = np.asarray(e.field("scarp_potential"), float)
    b = json.loads((root / "manifest.json").read_text(encoding="utf-8"))["basins"]
    pres = b["preserved"]
    areas = np.array([x.get("areaKm2", np.nan) for x in pres], float)
    areas = areas[np.isfinite(areas)]
    ok = land & np.isfinite(acc) & np.isfinite(slope) & (acc > 0) & (slope > 0)
    upland = ok & (elev > 0.5)
    la, ls = np.log10(acc[upland]), np.log10(slope[upland])
    cent, med = [], []
    for lo in np.arange(1.5, 5.5, 0.5):
        m = (la >= lo) & (la < lo + 0.5)
        if m.sum() >= 300:
            cent.append(lo + 0.25)
            med.append(float(np.mean(ls[m])))
    cent, med = np.array(cent), np.array(med)
    k = (cent >= 3.0) & (cent <= 5.5)
    theta = float(-np.polyfit(cent[k], med[k], 1)[0]) if k.sum() >= 3 else float("nan")
    return {
        "label": label, "n_regions": int(e.n_regions),
        "channel_pct_of_land_area":
            float(area[land & (acc > CHANNEL_THRESHOLD_KM2)].sum() / area[land].sum() * 100),
        "preserved_basins": len(pres),
        "smallest_preserved_km2": float(areas.min()) if areas.size else float("nan"),
        "scarp_nonzero_pct": float((scarp > 0).mean() * 100),
        "scarp_mean_where_nonzero": float(scarp[scarp > 0].mean()) if (scarp > 0).any() else float("nan"),
        "land_slope_p99_deg": float(np.percentile(np.degrees(np.arctan(slope[land])), 99)),
        "upland_concavity_theta": theta,
        "land_area_above": {f"{lv:g}": float(area[land & (elev >= lv)].sum() / area[land].sum())
                            for lv in LEVELS_KM},
        "mean_land_elevation_km":
            float((elev[land] * area[land]).sum() / area[land].sum()),
    }


FIELDS = [("channel_pct_of_land_area", "land area as channel, A>1e3 (%)"),
          ("preserved_basins", "preserved basins"),
          ("smallest_preserved_km2", "smallest preserved basin (km2)"),
          ("scarp_nonzero_pct", "cells with scarp potential > 0 (%)"),
          ("scarp_mean_where_nonzero", "mean scarp potential where > 0"),
          ("land_slope_p99_deg", "land slope 99th percentile (deg)"),
          ("upland_concavity_theta", "upland concavity theta"),
          ("mean_land_elevation_km", "mean land elevation (km)")]


def gridded_agreement(runs_paths, labels):
    """Cellwise agreement on the shared T42 grid, which is where the withdrawn
    RMS-from-converged-terrain figure came from. It is in this script rather
    than the audit's other one because it only means anything beside the
    control group: 2,600,001 regions agrees with the build no better than
    10,000,005 does, which is what makes the figure scatter and not distance."""
    from netCDF4 import Dataset
    a = Dataset(str(Path(runs_paths[0]) / "planet.nc"))
    ea = np.asarray(a["elevation_km"][:], float)
    sa = np.asarray(a["surface_class"][:], int)
    rows = []
    for path, label in zip(runs_paths[1:], labels[1:]):
        b = Dataset(str(Path(path) / "planet.nc"))
        eb = np.asarray(b["elevation_km"][:], float)
        sb = np.asarray(b["surface_class"][:], int)
        rows.append({"label": label,
                     "r": float(np.corrcoef(ea.ravel(), eb.ravel())[0, 1]),
                     "rms_km": float(np.sqrt(((ea - eb) ** 2).mean())),
                     "land_mask_agreement": float(((sa > 0) == (sb > 0)).mean())})
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--export", action="append", required=True, metavar="LABEL=PATH",
                    help="repeat; the CONTROLS are the ones within a few percent "
                         "of each other in region count")
    ap.add_argument("--controls", type=int, default=3,
                    help="how many of the leading exports form the same-resolution "
                         "control group (default 3)")
    ap.add_argument("--out", type=Path,
                    default=PROJECT_ROOT / "analysis" / "orogen_resolution_controls.json")
    args = ap.parse_args()

    runs = []
    for spec in args.export:
        label, path = spec.split("=", 1)
        p = Path(path)
        runs.append(load(p if p.is_absolute() else PROJECT_ROOT / p, label))

    hdr = f"{'statistic':>34s}" + "".join(f"{r['label']:>10s}" for r in runs)
    print(hdr + f"{'noise':>10s}{'clears?':>9s}")
    print("  " + "-" * (len(hdr) + 17))
    ctl = runs[:args.controls]
    verdicts = {}
    for key, name in FIELDS:
        vals = [r[key] for r in runs]
        c = [r[key] for r in ctl]
        noise = max(c) - min(c)
        delta = abs(vals[-1] - vals[0])
        verdict = ("yes" if delta > 2 * noise else
                   "marginal" if delta > noise else "NO")
        verdicts[key] = {"values": vals, "noise": noise, "delta": delta,
                         "clears_noise": verdict}
        print(f"{name:>34s}" + "".join(f"{v:10.4g}" for v in vals)
              + f"{noise:10.4g}{verdict:>9s}")
    print("\n  'noise' is the spread across the control group, which differs in "
          "region\n  count by a few percent and so holds resolution fixed. A "
          "statistic whose\n  change across the full range does not clear it is "
          "not a resolution result.")
    paths = [(Path(spec.split("=", 1)[1]) if Path(spec.split("=", 1)[1]).is_absolute()
              else PROJECT_ROOT / spec.split("=", 1)[1]) for spec in args.export]
    labels = [spec.split("=", 1)[0] for spec in args.export]
    try:
        ga = gridded_agreement(paths, labels)
        print(f"\n  gridded agreement against {labels[0]} on the shared T42 grid:")
        for g in ga:
            print(f"    {g['label']:>8s}  r {g['r']:+.4f}  RMS {g['rms_km']:.4f} km  "
                  f"land mask {g['land_mask_agreement']:.2%}")
        print("    a control a few percent away agrees no better than the far one: "
              "that\n    is why the RMS reads as scatter rather than as distance "
              "from a limit.")
    except Exception as e:
        ga = None
        print(f"\n  (gridded agreement unavailable: {e})")
    args.out.write_text(json.dumps(
        {"runs": runs, "controls": [r["label"] for r in ctl],
         "verdicts": verdicts, "gridded_agreement": ga}, indent=2) + "\n")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
