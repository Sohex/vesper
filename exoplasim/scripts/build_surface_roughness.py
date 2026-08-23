#!/usr/bin/env python3
"""Aerodynamic roughness length per cell, surface code 0173.

`configure()` clears every surface field when handed a landmap, so `dz0clim`
falls back to the uniform namelist default of `dz0land = 2.0 m` over all land.
That is the total roughness, and it is used directly: `NVEG = 0` in these runs,
so `landmod.f90:409` takes `dz0 = dz0clim` and SIMBA's separate vegetation and
orographic terms never run.

A single 2.0 m over all land is wrong in a way that matters here. It asserts
forest-scale roughness over the 16.5% of land that is salt crust and playa, where
real values are nearer 0.001 m, and those surfaces are not scattered: they are the
floors of closed basins, which are flat by construction and are exactly the cells
the carve verdict integrates evaporation over. Roughness enters the turbulent
exchange coefficient as

    ce = k^2 / ln(z_ref / z0)^2

so at this planet's 141.6 m lowest level the default gives 0.00882 against 0.00114
for a real playa surface: **7.7x too much exchange** over the most barren ground
on the planet, and therefore too much evaporation from precisely the basins whose
water balance decides whether they survive.

## What this computes

Total roughness combines a surface term and a subgrid orographic term in
quadrature, which is the same combination `simba.f90:474` uses when it builds
`dz0` from vegetation and topography:

    z0 = sqrt(z0_surface^2 + z0_orographic^2)

`z0_surface` is per mesh region, from land cover: barren classes take a bare-ground
value, everything else takes a canopy value scaled by forest fraction. It is
integrated to the grid over land only, so a coastal cell is not dragged toward
open water.

`z0_orographic` is the part worth retaining the native mesh for. The 10M
fine-support reference contributes about 1,221 regions per global T42 cell on
average (and about 76 even at T170), enough to measure the distribution of
elevation *inside* the cell rather than infer it from resolved slope. Counts
over land vary and are reported by the builder; SPAT-3 is the shared artifact
for comparisons across the full T21/T42/T85/T127/T170 ladder. That subgrid
relief is what a gridded mean elevation cannot supply.

## The one free parameter, and why it is anchored rather than chosen

The constant relating subgrid relief to an effective roughness is not something
this project can derive; formulations in the literature differ by more than an
order of magnitude. Rather than pick one, it is solved so that the **area-weighted
land mean of the result equals the model's own `dz0land`**. That keeps the global
value ExoPlaSim was tuned against and changes only the distribution, which is the
part that is physically wrong. `--target-mean` overrides it if there is ever a
reason to move the mean as well, and the solved coefficient is reported so the
implied relation can be judged.

This is deliberately the same move `build_surface_albedo.py --mode scaled` makes:
where a global constant is calibrated and its spatial pattern is not, keep the
constant and supply the pattern.

**This changes climate results.** It is off unless `model.roughness_source`
is set in `config/planet.yaml` (currently `lithology`); changing the key is a
configuration change `continue_exoplasim.py` compares for and refuses to
resume across, so moving it blocks any run in flight.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import yaml

from _paths import CONFIG, INPUTS, PROJECT_ROOT
from sra import write_sra
from builds import resolution_of, grid_export, mesh_export
# ONE SOURCE for how much forest the land carries. `model.land_albedo_source`
# already decides it, and code 212 is derived from that decision, so importing
# the mapping is what keeps this field and the albedo field describing one land
# cover. Re-deriving it here would be two formulations of one quantity, which
# is failure class 17. The same idiom continue_exoplasim.py uses on
# run_exoplasim.py, and for the same reason.
from build_surface_albedo import MODE_FOREST_FRACTION
from gridding import land_weighted, region_cells
from provenance import config_stamp
from orogen import Export, LAND

ROUGHNESS_CODE = 173

# landmod.f90:51, the uniform default this replaces.
EXOPLASIM_DZ0LAND_M = 2.0

# Ocean roughness is set by the model from wind, so the value written on ocean
# cells is inert. It is written as open water's own scale rather than zero so the
# field is never mistaken for a land-only mask.
OCEAN_Z0_M = 1.5e-4

# Surface roughness by cover, metres. Bare ground spans 1e-4 (salt crust) to
# 1e-2 (rocky desert); 3e-3 sits in that range and well below anything vegetated.
DEFAULT_BARE_Z0_M = 0.003
DEFAULT_CANOPY_Z0_M = 0.25       # closed grass to open woodland
DEFAULT_FOREST_Z0_M = 1.0        # closed forest


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--mesh", type=Path, default=None)
    ap.add_argument("--grid", type=Path, default=None)
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--target-mean", type=float, default=None,
                    help="area-weighted land mean to anchor to; defaults to "
                         "ExoPlaSim's own dz0land")
    ap.add_argument("--bare-z0", type=float, default=DEFAULT_BARE_Z0_M)
    ap.add_argument("--canopy-z0", type=float, default=DEFAULT_CANOPY_Z0_M)
    ap.add_argument("--forest-z0", type=float, default=DEFAULT_FOREST_Z0_M)
    ap.add_argument("--forest-fraction", type=float, default=None,
                    help="override the fraction implied by "
                         "model.land_albedo_source, the same override "
                         "build_surface_albedo.py takes, for a pair that "
                         "moves both fields together")
    ap.add_argument("--lakes", type=Path, default=None,
                    help="surface_water.nc; open water is smooth, so lake "
                         "regions take the ocean value before integration")
    args = ap.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    model = config["model"]
    # Deliberately from the grid, not from config: see builds.resolution_of.
    mesh = Export(args.mesh or mesh_export(config))
    grid_dir = args.grid or grid_export(config)
    # From the grid, not from config: the two differ exactly when someone
    # builds for another resolution, which is when the filename matters.
    resolution = resolution_of(grid_dir)

    rock = mesh.field("substrate_class").astype(int)
    is_land = mesh.surface_class == LAND
    area = mesh.cell_area.astype(np.float64)
    elev_km = mesh.field("elevation_km").astype(np.float64)

    codes = [r["code"] for r in mesh.manifest["lithology"]["rockClasses"]]
    barren = np.zeros(rock.shape, dtype=bool)
    for code in model.get("barren_rock_classes", []):
        if code in codes:
            barren |= rock == codes.index(code)

    # Surface term. The forest fraction comes from `model.land_albedo_source`,
    # the same key build_surface_albedo.py resolves its mode from, so roughness
    # and code 212 describe ONE land cover. It used to read
    # `model.forest_fraction_assumed`, a key the config does not carry, so it
    # silently used zero forest while 212 asserted half a canopy; that is
    # CLIM-36. An explicit --forest-fraction still overrides, for a sensitivity
    # pair. `modelled` is absent from the mapping because it does not imply a
    # fraction, it reads tree cover per cell, so it falls back to no blend here
    # until this reads that field too.
    mode = str(model.get("land_albedo_source", "lithology"))
    forest_fraction = (args.forest_fraction if args.forest_fraction is not None
                       else MODE_FOREST_FRACTION.get(mode))
    forest_fraction = float(forest_fraction) if forest_fraction else None
    canopy_z0 = args.canopy_z0
    if forest_fraction is not None:
        canopy_z0 = ((1.0 - forest_fraction) * args.canopy_z0
                     + forest_fraction * args.forest_z0)
    z0_surface = np.where(barren, args.bare_z0, canopy_z0)

    if args.lakes is not None:
        from netCDF4 import Dataset
        with Dataset(args.lakes) as lds:
            lake = np.asarray(lds["lake"][:]).astype(bool)
            lake_terrain = getattr(lds, "terrain_hash", None)
        if lake_terrain and lake_terrain != mesh.terrain_hash:
            raise SystemExit(
                f"lake solution is on terrain {lake_terrain[:16]}, mesh is "
                f"{mesh.terrain_hash[:16]}; re-run surface_water.py")
        z0_surface = np.where(lake & is_land, OCEAN_Z0_M, z0_surface)

    _, z0_surf_grid, _ = land_weighted(mesh, grid_dir, z0_surface)

    # Orographic term: the standard deviation of elevation among the mesh regions
    # inside each cell, which is subgrid relief by construction.
    cells, nlat, nlon = region_cells(mesh, grid_dir)
    n = nlat * nlon
    sel = is_land
    w = np.bincount(cells[sel], weights=area[sel], minlength=n)
    m1 = np.bincount(cells[sel], weights=area[sel] * elev_km[sel], minlength=n)
    m2 = np.bincount(cells[sel], weights=area[sel] * elev_km[sel] ** 2, minlength=n)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(w > 0, m1 / np.maximum(w, 1e-30), 0.0)
        var = np.where(w > 0, m2 / np.maximum(w, 1e-30) - mean ** 2, 0.0)
    sigma_m = np.sqrt(np.maximum(var, 0.0)).reshape(nlat, nlon) * 1000.0
    counts = np.bincount(cells[sel], minlength=n).reshape(nlat, nlon)

    land_cells = w.reshape(nlat, nlon) > 0
    gw = np.fromfile(grid_dir / "grid" / "gauss_weights.bin", dtype="float64")
    weight = np.where(land_cells, w.reshape(nlat, nlon), 0.0)

    def land_mean(field):
        return float((field * weight).sum() / weight.sum())

    # Solve the orographic coefficient so the land mean lands on the target.
    target = args.target_mean if args.target_mean is not None else EXOPLASIM_DZ0LAND_M
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        trial = land_mean(np.sqrt(z0_surf_grid ** 2 + (mid * sigma_m) ** 2))
        if trial < target:
            lo = mid
        else:
            hi = mid
    k_oro = 0.5 * (lo + hi)
    z0 = np.sqrt(z0_surf_grid ** 2 + (k_oro * sigma_m) ** 2)
    field = np.where(land_cells, z0, OCEAN_Z0_M)

    output = args.output or (INPUTS / resolution.lower()
                             / f"orogen_{resolution}_surf_{ROUGHNESS_CODE:04d}.sra")
    output.parent.mkdir(parents=True, exist_ok=True)
    write_sra(output, ROUGHNESS_CODE, field)

    ce = lambda z0v: 0.16 / np.log(141.6 / np.maximum(z0v, 1e-6)) ** 2
    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "code": ROUGHNESS_CODE,
        "field": "dz0clim, total aerodynamic roughness length, metres",
        "terrain_hash": mesh.terrain_hash,
        "replaces_uniform": EXOPLASIM_DZ0LAND_M,
        "surface_z0_m": {"bare": args.bare_z0, "canopy": canopy_z0,
                         "forest": args.forest_z0,
                         "forest_fraction_used": forest_fraction},
        "orographic": {
            "coefficient_solved": round(k_oro, 6),
            "relation": "z0_oro = coefficient * stdev(elevation) within the cell",
            "subgrid_stdev_m": {
                "min": round(float(sigma_m[land_cells].min()), 2),
                "median": round(float(np.median(sigma_m[land_cells])), 2),
                "max": round(float(sigma_m[land_cells].max()), 2)},
            "mesh_regions_per_land_cell_median": int(np.median(counts[land_cells])),
        },
        "land_mean_m": round(land_mean(z0), 5),
        "land_min_m": round(float(z0[land_cells].min()), 6),
        "land_max_m": round(float(z0[land_cells].max()), 4),
        "exchange_coefficient": {
            "note": "ce = k^2 / ln(141.6 / z0)^2, the quantity roughness acts through",
            "uniform_default": round(float(ce(EXOPLASIM_DZ0LAND_M)), 6),
            # ce rises with z0, so the roughest cell carries the largest ce.
            "this_field_land_max": round(float(ce(z0[land_cells].max())), 6),
            "this_field_land_min": round(float(ce(z0[land_cells].min())), 6),
        },
        "file": str(output),
    }
    report_path = output.parent / f"roughness_{resolution.lower()}_report.json"
    # Provenance stamp; lib/provenance.py owns the shape and the inert set.
    report.update(config_stamp(config, "exoplasim/scripts/build_surface_roughness.py"))
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"land-mean z0 {land_mean(z0):.4f} m against the uniform "
          f"{EXOPLASIM_DZ0LAND_M} m it replaces")
    print(f"  range {z0[land_cells].min():.5f} to {z0[land_cells].max():.3f} m")
    print(f"  orographic coefficient solved to {k_oro:.5f}, subgrid stdev median "
          f"{np.median(sigma_m[land_cells]):.1f} m")
    print(f"  exchange coefficient {ce(z0[land_cells].min()):.5f} (flattest, "
          f"barest) to {ce(z0[land_cells].max()):.5f} (roughest), against a "
          f"uniform {ce(EXOPLASIM_DZ0LAND_M):.5f}")
    print(f"wrote {output}\n      {report_path}")


if __name__ == "__main__":
    main()
