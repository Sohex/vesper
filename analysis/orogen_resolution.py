#!/usr/bin/env python3
"""Is Orogen's terrain converged in the mesh, and what does a finer one buy?

    python analysis/orogen_resolution.py --export "2.5M build=source/precarve-craton/exoplasim-T42" \
                                         --export "10M=/scratch/n10m"

Worldbuilding. Vesper is an invented planet and this script measures the
generator that makes its terrain: whether the elevation field a World Orogen
export carries is a property of the world or of the mesh it was sampled on.
Every quantity here is a modelled field of an invented planet.

`notes/audits/orogen-resolution.md` is the argument; this is the measurement
made re-runnable. It writes `analysis/orogen_resolution.json`.

## Three questions, and they do not have the same answer

RELIEF asks how much elevation varies between two points a given distance
apart. The obvious estimator, the standard deviation of elevation inside a cap
of radius r, punishes a coarse mesh twice: it holds fewer cells per cap, and a
standard deviation from three samples is biased low, so the coarse mesh scores
low partly for having less relief and partly for having less data. Pooling
every pair removes that:

    gamma(r) = 0.5 * mean[(z_i - z_j)^2] over land pairs separated by at most r

with sqrt(gamma) in the units of the field. This is what settles.

HYPSOMETRY asks how much land AREA sits above a level. It is a different
statistic from relief, it is area-weighted rather than per-cell, and it is the
one that does not settle: a coarse cell gives one averaged value to a large
area, so ground high in only part of that cell is recorded as high throughout
it. Reporting relief alone would call the terrain converged and be wrong.

SANITY asks whether the export is well formed at all, since the build already
sits at the generator's design ceiling of 2,560,000 regions and anything above
it is outside the range the tool was tuned in. The closure identity is
`cell_area` summing to the sphere: the centroidal dual overshoots by a fixed
1.00068 and that number is the check, not 1.0. See
`notes/audits/mesh-dual-area.md`.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "lib"))

from lib.orogen import Export  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

LAG_KM = [5.0, 10.0, 15.19, 25.0, 50.0, 100.0, 200.0]
LEVELS_KM = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
N_CENTRES = 6000
SEED = 20260820


def load(root: Path) -> dict:
    e = Export(root, require_known_build=False)
    xyz = np.column_stack([e.x, e.y, e.z]).astype(np.float64)
    xyz /= np.linalg.norm(xyz, axis=1, keepdims=True)
    out = {
        "export": e,
        "n_regions": int(e.n_regions),
        "radius_km": float(e.radius_km),
        "xyz": xyz,
        "elev_km": np.asarray(e.elevation_km, float),
        "area_km2": np.asarray(e.cell_area, float),
        "land": np.asarray(e.surface_class, int) > 0,
    }
    try:
        out["pre"] = np.asarray(e.field("elevation_pre_erosion"), float)
    except Exception:
        out["pre"] = None
    return out


def sanity(d: dict) -> dict:
    R, n = d["radius_km"], d["n_regions"]
    sphere = 4 * np.pi * R * R
    area_sum = float(d["area_km2"].sum())
    el, ar, land = d["elev_km"], d["area_km2"], d["land"]
    return {
        "n_regions": n,
        "mean_edge_km": float(np.pi * R / np.sqrt(n)),  # Orogen avgEdgeKm; the project's spacing number
        "cell_area_sum_km2": area_sum,
        "cell_area_over_sphere": area_sum / sphere,
        "non_finite_elevations": int((~np.isfinite(el)).sum()),
        "elevation_min_km": float(el.min()),
        "elevation_max_km": float(el.max()),
        "land_fraction_by_count": float(land.mean()),
        "land_fraction_by_area": float(ar[land].sum() / area_sum),
        "unit_sphere_deviation": float(np.abs(np.linalg.norm(d["xyz"], axis=1) - 1).max()),
    }


def hypsometry(d: dict) -> dict:
    el, ar, land = d["elev_km"], d["area_km2"], d["land"]
    e, a = el[land], ar[land]
    tot = float(a.sum())
    return {
        "land_area_km2": tot,
        "area_weighted_mean_land_elevation_km": float((e * a).sum() / tot),
        "area_fraction_above": {f"{lv:g}": float(a[e >= lv].sum() / tot) for lv in LEVELS_KM},
    }


def variogram(d: dict, rng: np.random.Generator) -> dict:
    """sqrt of the semivariance at each lag, pooled over all land pairs."""
    R = d["radius_km"]
    tree = cKDTree(d["xyz"])
    land_idx = np.flatnonzero(d["land"])
    centres = rng.choice(land_idx, size=min(N_CENTRES, land_idx.size), replace=False)
    fields = {"elevation_km": d["elev_km"]}
    if d["pre"] is not None:
        fields["elevation_pre_erosion"] = d["pre"]
    out: dict = {}
    for lag in LAG_KM:
        chord = 2 * np.sin(lag / R / 2.0)
        neighbours = tree.query_ball_point(d["xyz"][centres], chord, workers=-1)
        per_field = {}
        for name, a in fields.items():
            ssq, npair = 0.0, 0
            for c, nb in zip(centres, neighbours):
                if len(nb) < 2:
                    continue
                diff = a[nb] - a[c]
                ssq += float(diff @ diff)
                npair += len(nb) - 1
            per_field[name] = {
                "sqrt_semivariance": float(np.sqrt(0.5 * ssq / npair)) if npair else None,
                "pairs": npair,
            }
        out[f"{lag:g}"] = per_field
    return out


def _sha256(path: Path) -> str | None:
    if not path.is_file():
        return None
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def provenance(exports: dict[str, Path]) -> dict:
    try:
        rev = subprocess.run(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT,
                             capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        rev = None
    return {
        "generated": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "git_rev": rev,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "script_sha256": _sha256(Path(__file__).resolve()),
        "inputs": {
            label: {"path": str(root),
                    "manifest_sha256": _sha256(root / "manifest.json")}
            for label, root in exports.items()
        },
        "settings": {"lags_km": LAG_KM, "levels_km": LEVELS_KM,
                     "n_centres": N_CENTRES, "seed": SEED},
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--export", action="append", required=True, metavar="LABEL=PATH",
                    help="an Orogen export to measure; repeat to compare region counts")
    ap.add_argument("--out", type=Path,
                    default=PROJECT_ROOT / "analysis" / "orogen_resolution.json")
    args = ap.parse_args()

    exports: dict[str, Path] = {}
    for spec in args.export:
        if "=" not in spec:
            ap.error(f"--export wants LABEL=PATH, got {spec!r}")
        label, path = spec.split("=", 1)
        p = Path(path)
        exports[label] = p if p.is_absolute() else (PROJECT_ROOT / p)

    results: dict = {}
    for label, root in exports.items():
        if not root.exists():
            ap.error(f"{label}: {root} does not exist")
        d = load(root)
        results[label] = {
            "sanity": sanity(d),
            "hypsometry": hypsometry(d),
            "variogram": variogram(d, np.random.default_rng(SEED)),
        }
        s, h = results[label]["sanity"], results[label]["hypsometry"]
        print(f"=== {label}: {s['n_regions']:,} regions, "
              f"{s['mean_edge_km']:.2f} km mean edge ===")
        print(f"  cell_area / sphere      {s['cell_area_over_sphere']:.10f}  "
              f"(centroidal dual overshoots; 1.0 is NOT the expected value)")
        print(f"  non-finite elevations   {s['non_finite_elevations']}")
        print(f"  land fraction by area   {s['land_fraction_by_area']:.4f}")
        print(f"  mean land elevation     {h['area_weighted_mean_land_elevation_km']:.4f} km")

    labels = list(results)
    print(f"\nsqrt(semivariance) of elevation_km vs lag, land pairs")
    print(f"  {'lag km':>8s}" + "".join(f"{l:>16s}" for l in labels))
    for lag in LAG_KM:
        row = f"  {lag:8.2f}"
        for l in labels:
            v = results[l]["variogram"][f"{lag:g}"]["elevation_km"]["sqrt_semivariance"]
            row += f"{v:16.4f}" if v is not None else f"{'--':>16s}"
        print(row)

    print(f"\nfraction of land AREA above each level")
    print(f"  {'level km':>8s}" + "".join(f"{l:>16s}" for l in labels))
    for lv in LEVELS_KM:
        row = f"  {lv:8.1f}"
        for l in labels:
            row += f"{results[l]['hypsometry']['area_fraction_above'][f'{lv:g}']:16.5f}"
        print(row)

    payload = {"provenance": provenance(exports), "results": results}
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
