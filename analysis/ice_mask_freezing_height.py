#!/usr/bin/env python3
"""Where ice can persist on this planet, at the mesh's resolution rather than the model's.

    python analysis/ice_mask_freezing_height.py \
        --climatology exoplasim/analysis/climatology/bootstrap_regular_climatology.nc

Worldbuilding. Vesper is an invented planet; this is about the ice in a toy
climate model of it and the terrain a separate generator makes.

PHYS-13's evidence, made re-runnable. Orogen places the ice that carves this
terrain by latitude and dimensionless elevation alone, with Earth-tuned
constants, so at the build's `glacialErosion` of 0.8 the ice line sits at 58
degrees put there by the erosion slider and nothing else -- no obliquity, no
orbit, no stellar spectrum.

The obvious fix, driving it from the climatology's `glac`, does not work: that
field is identically ZERO over every cell and every output bin while the run
manifest records `"glaciers": true`. The model ran and formed no permanent land
ice, and that is a fact about the GRID rather than the planet. A cell's
temperature is evaluated at its MEAN elevation, and the mean of a coarse cell is
far below the ground inside it, so the model never sees the cold ground.

What works is to take the smooth quantity the grid genuinely resolves, the
warmest-month surface temperature, and correct it to each mesh region's own
elevation with the measured lapse rate. `maps/build_basemap.py` already does
exactly this for the basemap's ice shading and `lib/lapse.py` names it the
freezing-height criterion.

The criterion is THERMAL: it says where ice can persist, not where a glacier
forms, which additionally needs accumulation.

`--write-mask PATH` emits the mask Orogen consumes: one float32 per mesh region
in region order, plus a `PATH.json` sidecar. The sidecar is not optional. The
mask is matched to the terrain BY REGION INDEX, and a mesh is a pure function of
the seed and the region count, so those two are the identity the generator
checks before indexing anything. Matching by longitude across this boundary has
silently matched zero cells three times; see CLAUDE.md rule 3.

A mask carries no duration. It says WHERE the model's climate keeps ground below
freezing in the warmest bin, never for how long, so the glacial strength slider
Orogen scales it by stays a declared choice; `docs/src/reference/no-time-axis.md`.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "lib"))

from lib.orogen import Export  # noqa: E402
import gridding  # noqa: E402
import lapse  # noqa: E402
from lib import builds  # noqa: E402

FREEZE_K = 273.15


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--climatology", type=Path, required=True)
    ap.add_argument("--build", type=Path, default=None,
                    help="export root carrying raw/, the per-region mesh arrays; "
                         "defaults to the configured source_build")
    ap.add_argument("--grid-export", type=Path, default=None,
                    help="export whose GRID the climatology was run on. Only the T42 "
                         "export keeps raw/, so a climatology at another truncation "
                         "needs the mesh from one export and the grid from another; "
                         "the region-to-cell mapping is built from this one")
    ap.add_argument("--margins-k", type=float, nargs="*", default=[0.0, -5.0, -10.0],
                    help="report the area below freezing plus each margin, as a "
                         "crude sensitivity to the criterion's sharpness")
    ap.add_argument("--out", type=Path,
                    default=PROJECT_ROOT / "analysis" / "ice_mask_freezing_height.json")
    ap.add_argument("--write-mask", type=Path, default=None,
                    help="also write the per-region mask Orogen's --ice-mask reads, "
                         "plus its .json sidecar")
    ap.add_argument("--mask-margin-k", type=float, default=0.0,
                    help="offset on the freezing criterion used for the WRITTEN mask, "
                         "declared before the mask is made rather than tuned after")
    a = ap.parse_args()

    root = a.build or builds.mesh_export()
    grid_root = Path(a.grid_export) if a.grid_export else Path(root)
    e = Export(Path(root), require_known_build=False)
    d = Dataset(str(a.climatology))
    ts = np.asarray(d["ts"][:], float)
    warmest = ts.max(axis=0)
    clim_lat = np.asarray(d["lat"][:], float)
    rate = lapse.environmental_lapse_k_per_km(climatology_file=a.climatology,
                                              season="warmest")

    # rule 3: the export and the climatology label their columns differently, so
    # this goes by INDEX through the one grid convention and never by longitude.
    row, col = gridding.climatology_cells(e, grid_root, clim_lat)
    land = np.asarray(e.surface_class, int) > 0
    area = np.asarray(e.cell_area, float)
    elev = np.asarray(e.elevation_km, float)

    gnc = Dataset(str(grid_root / "planet.nc"))
    grid_elev = np.asarray(gnc["elevation_km"][:], float)
    t_cell = warmest[row, col]
    z_cell = grid_elev[row, col]
    t_local = t_cell - rate * (elev - z_cell)

    la = area[land].sum()
    glac = np.asarray(d["glac"][:], float) if "glac" in d.variables else None
    out = {
        "climatology": str(a.climatology),
        "mesh_export": str(root),
        "grid_export": str(grid_root),
        "lapse_k_per_km_warmest": rate,
        "grid_cells_with_glac": (int((glac.max(axis=0) > 0).sum())
                                 if glac is not None else None),
        "uncorrected_land_below_freezing":
            float(area[land & (t_cell < FREEZE_K)].sum() / la),
        "corrected_land_below_freezing":
            float(area[land & (t_local < FREEZE_K)].sum() / la),
        "by_margin_k": {},
    }
    # what the grid cannot see: the ground inside a cell against its mean
    if all(k in gnc.variables for k in ("orog_max", "orog_mean", "surface_class")):
        gm = np.asarray(gnc["orog_mean"][:], float)
        gx = np.asarray(gnc["orog_max"][:], float)
        gl = np.asarray(gnc["surface_class"][:], int) > 0
        exc = (gx - gm)[gl]
        out["subgrid_peak_excess_km"] = {
            "mean": float(exc.mean()), "median": float(np.median(exc)),
            "p90": float(np.percentile(exc, 90)), "max": float(exc.max())}
        out["temperature_the_grid_never_sees_k"] = {
            "mean": float(exc.mean() * rate), "p90": float(np.percentile(exc, 90) * rate)}

    print(f"warm-season environmental lapse rate: {rate:.3f} K/km")
    if glac is not None:
        print(f"grid cells with any glac: {out['grid_cells_with_glac']} "
              f"of {glac.shape[1] * glac.shape[2]}")
    if "subgrid_peak_excess_km" in out:
        s = out["subgrid_peak_excess_km"]
        print(f"sub-grid peak excess on land: mean {s['mean']:.3f} km, "
              f"p90 {s['p90']:.3f} km")
        t = out["temperature_the_grid_never_sees_k"]
        print(f"  which is {t['mean']:.1f} K the grid runs too warm on its own high "
              f"ground, {t['p90']:.1f} K at p90")
    print(f"\nland below freezing in the warmest month:")
    print(f"  on the model's own orography   {out['uncorrected_land_below_freezing']:7.3%}")
    print(f"  corrected to the mesh          {out['corrected_land_below_freezing']:7.3%}")
    for m in a.margins_k:
        v = float(area[land & (t_local < FREEZE_K + m)].sum() / la)
        out["by_margin_k"][f"{m:g}"] = v
        print(f"     margin {m:+5.1f} K              {v:7.3%}")
    if a.write_mask is not None:
        write_mask(a.write_mask, root, e, t_local, land, a.mask_margin_k, rate, a.climatology)
        out["mask"] = {"path": str(a.write_mask), "margin_k": a.mask_margin_k}

    a.out.write_text(json.dumps(out, indent=2) + "\n")
    print(f"\nwrote {a.out}")


def write_mask(path: Path, root: Path, e, t_local, land, margin_k: float,
               rate: float, climatology: Path) -> None:
    """Emit the per-region ice mask and the sidecar that identifies its mesh.

    The mask is BINARY at the criterion rather than a ramp across it. A ramp
    would need a width, and no measurement here sets one; a declared width
    chosen after seeing the map is not a criterion. The sensitivity to the
    criterion's sharpness is reported instead, by `--margins-k`.

    Ocean regions are zero. Orogen zeroes them again on its own ocean mask, so
    this is agreement rather than reliance: the two disagree over dry
    closed-basin floor below sea level, and land here is `surface_class`, which
    is the side that keeps it (CLAUDE.md rule 1).
    """
    manifest = json.loads((Path(root) / "manifest.json").read_text())
    values = np.zeros(t_local.shape, dtype=np.float32)
    values[land & (t_local < FREEZE_K + margin_k)] = 1.0

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(values.tobytes())
    sidecar = {
        "numRegions": int(values.size),
        "seed": manifest.get("seed"),
        "regionsRequested": manifest.get("params", {}).get("N"),
        "sourceBuild": str(root),
        "terrainHash": manifest.get("hashes", {}).get("finalElevation"),
        "climatology": str(climatology),
        "criterion": "warmest-bin surface temperature, lapse-corrected to the "
                     "mesh region's own elevation, below freezing",
        "marginK": margin_k,
        "lapseKPerKmWarmest": rate,
        "glaciatedRegions": int(values.sum()),
        "generator": "analysis/ice_mask_freezing_height.py",
    }
    Path(str(path) + ".json").write_text(json.dumps(sidecar, indent=2) + "\n")
    print(f"\nwrote {path} ({sidecar['glaciatedRegions']} glaciated regions "
          f"of {sidecar['numRegions']}) and its sidecar")


if __name__ == "__main__":
    main()
