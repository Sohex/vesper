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

with `k` von Karman's constant and `z_ref` the height of the lowest model level,
which `lib/lapse.py:reference_height_m` derives hypsometrically at THIS planet's
gravity. Over the liquid-water span the default gives about 8.1 to 8.5 times the
exchange a real playa surface would, and therefore that much too much evaporation
from precisely the basins whose water balance decides whether they survive. The
contrast is reported bracketed rather than at one temperature because `z_ref` is
linear in the air temperature and this step runs before any climatology exists to
measure one from.

## What this computes

Total roughness combines a surface term and a subgrid orographic term in
quadrature, which is the same combination `simba.f90:474` uses when it builds
`dz0` from vegetation and topography:

    z0 = sqrt(z0_surface^2 + z0_orographic^2)

`z0_surface` is per mesh region, from land cover: barren classes take a bare-ground
value, everything else takes a canopy value scaled by forest fraction. It is
integrated to the grid over land only, so a coastal cell is not dragged toward
open water.

**It is integrated in `ce`, not in the length.** The model applies one exchange
coefficient to a whole cell and the turbulent flux is linear in it, so what the
cell owes the atmosphere is the area mean of `ce` over its own surfaces, and the
roughness written out is the length that reproduces that mean. `ce` is
logarithmic in `z0`, so the two orders are different reductions, and this
world's land is the case where they part: bare ground and canopy are two orders
apart in `z0`, so a cell that is mostly playa with a canopy minority has its
exchange set by the minority once the lengths are mixed. SPAT-7 measured the
difference on the 10M mesh -- negligible in the land mean, 24 times inside the
`z_ref` bracket this file already declares, but past that bracket on a few
percent of land area which GROWS with refinement, and reaching 30% on closed
basin floors. Those cells are the reason this field exists. The inversion needs
`z_ref`, which is unknown before a climatology; it is taken at the midpoint of
the liquid-water bracket and the report carries what the whole bracket is worth,
which is 1.5e-4 of the land mean.

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

**A LADDER COMPARISON MUST PASS ONE COEFFICIENT TO EVERY RUNG.** The solve is per
grid, and it has to be, because the subgrid relief it multiplies is per grid:
measured on the 10M mesh the median falls by 5.96 from T21 to T170 while the
solved coefficient rises by 1.98 to hold the land mean. So two rungs built with
the defaults differ by their terrain AND by their calibration, and a convergence
claim taken across them is measuring both. `--orographic-coefficient` fixes it;
the default still solves, so a single-rung build is unchanged. This is SPAT-8's
constraint and `notes/audits/nonlinear-spatial-reductions.md` carries the
numbers.

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
from gridding import (cell_expectation, cell_moments, region_cells,
                      transfer_ledger)
from provenance import config_stamp
from orogen import Export, LAND
# ONE derivation of the height a bulk transfer coefficient is taken over.
# hydrography/scripts/carve_verdict.py builds the same z_ref for the same
# reason; this file used to carry the EARTH-gravity answer as a literal.
from lapse import reference_height_m

ROUGHNESS_CODE = 173

KARMAN = 0.4                     # von Karman's constant; `ce` uses its square.

# The reference height is linear in the lowest-level air temperature and this
# step runs before any climatology exists to measure one from, so the reported
# exchange coefficients are BRACKETED over the range in which surface water is
# liquid -- which is the range the field's consumer, lake and playa
# evaporation, is defined over. Fixed here, before any field is built.
CE_BRACKET_K = (273.15, 313.15)

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
    ap.add_argument("--orographic-coefficient", type=float, default=None,
                    help="fix the coefficient relating subgrid relief to a "
                         "roughness instead of solving it on this grid. A "
                         "ladder comparison MUST pass one value to every rung: "
                         "the solve is per grid, so two rungs otherwise differ "
                         "by their terrain and by their calibration at once")
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

    # Orographic term: the standard deviation of elevation among the mesh regions
    # inside each cell, which is subgrid relief by construction. `cell_moments`
    # in `lib/gridding.py` is the operator; the three bincounts used to be here.
    cells, nlat, nlon = region_cells(mesh, grid_dir)
    n = nlat * nlon
    sel = is_land
    _, var, counts, land_flat = cell_moments(cells, n, area, elev_km, sel)
    w = np.bincount(cells[sel], weights=area[sel], minlength=n)
    sigma_m = np.sqrt(var).reshape(nlat, nlon) * 1000.0
    counts = counts.reshape(nlat, nlon)

    land_cells = land_flat.reshape(nlat, nlon)
    weight = np.where(land_cells, w.reshape(nlat, nlon), 0.0)

    def land_mean(field):
        return float((field * weight).sum() / weight.sum())

    # THE SURFACE TERM IS REDUCED IN THE EXCHANGE COEFFICIENT, NOT IN THE LENGTH.
    #
    # The model applies one `ce` to the whole cell and the turbulent flux is
    # linear in it, so what a cell owes the atmosphere is the AREA MEAN OF `ce`
    # over its own surfaces. `ce` is not linear in `z0` -- it is
    # `k^2 / ln(z_ref/z0)^2` -- so averaging the lengths first and taking `ce`
    # of that is not the same reduction, and this world's land is exactly the
    # case where the two part company: the barren classes carry a roughness two
    # orders below the canopy, so a cell that is mostly playa with a canopy
    # minority has its exchange set by the minority once the lengths are mixed.
    # SPAT-7 measured it on the 10M mesh: negligible in the land mean, 24 times
    # inside the step's own `z_ref` bracket, but past that bracket on a few
    # percent of land area which GROWS with refinement, and reaching 30% on the
    # closed-basin floors the carve verdict integrates evaporation over. Those
    # cells are the reason this field exists. `analysis/spatial_reduction_gap.py`
    # is the measurement and `notes/audits/nonlinear-spatial-reductions.md` the
    # finding.
    #
    # So: `ce` per mesh region, area-mean over the cell's land, and the
    # roughness written out is the length that reproduces that mean. The
    # inversion needs `z_ref`, which is not known before a climatology exists,
    # so it is taken at the MIDPOINT of the same liquid-water bracket the report
    # already declares and the bracket's own effect on the answer is reported
    # rather than hidden.
    z_ref_anchor = reference_height_m(0.5 * (CE_BRACKET_K[0] + CE_BRACKET_K[1]),
                                      config)

    def ce_of(z0v, z_ref: float):
        """Neutral bulk exchange coefficient at the lowest model level."""
        return KARMAN ** 2 / np.log(z_ref / np.maximum(z0v, 1e-6)) ** 2

    def region_z0(k_oro: float) -> np.ndarray:
        """Total roughness PER MESH REGION, the two terms in quadrature."""
        return np.sqrt(z0_surface ** 2
                       + (k_oro * sigma_m.reshape(-1)[cells]) ** 2)

    def reduce_ce(k_oro: float, z_ref: float):
        """`(cell z0, ce expectation, ce of the mean length, covered)`.

        The operator is `gridding.cell_expectation`, the NONLINEAR one, and it
        is called rather than reimplemented: it takes the LAW as a callable, so
        a caller cannot hand it a field that has already been reduced. The
        three bincounts that used to be here were the same reduction written a
        second time, which is what `lib/gridding.py` exists to stop.
        """
        ce_bar, z0_bar, covered = cell_expectation(
            cells, n, area, lambda v: ce_of(v, z_ref), region_z0(k_oro), sel)
        out = np.zeros(n)
        np.multiply(z_ref, np.exp(-KARMAN / np.sqrt(np.maximum(ce_bar, 1e-30))),
                    out=out, where=covered)
        return out.reshape(nlat, nlon), ce_bar, ce_of(z0_bar, z_ref), covered

    def effective_z0(k_oro: float, z_ref: float) -> np.ndarray:
        """Cell roughness whose `ce` is the area mean of the regions' own."""
        return reduce_ce(k_oro, z_ref)[0]

    # Solve the orographic coefficient so the land mean lands on the target.
    # 60 bisections resolve the coefficient to 1e-18 on [0, 1]; the 200 this
    # carried were free when the trial was arithmetic on the grid and are not
    # now that each one reduces the mesh.
    target = args.target_mean if args.target_mean is not None else EXOPLASIM_DZ0LAND_M
    if args.orographic_coefficient is not None:
        k_oro = float(args.orographic_coefficient)
    else:
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if land_mean(effective_z0(mid, z_ref_anchor)) < target:
                lo = mid
            else:
                hi = mid
        k_oro = 0.5 * (lo + hi)
    z0, ce_expectation, ce_of_mean, covered = reduce_ce(k_oro, z_ref_anchor)
    field = np.where(land_cells, z0, OCEAN_Z0_M)

    # THE JENSEN GAP IS THE FIELD'S OWN EVIDENCE THAT THE ORDER MATTERS.
    # `expectation - law(mean)` in `ce`, relative to the expectation, per cell.
    # It is what separates this reduction from the one this file used to do,
    # and reporting it means a later reader does not have to re-derive SPAT-7's
    # measurement to see whether the two orders part on THIS build and THIS
    # grid. A gap that reads zero everywhere would say the correction is inert
    # here, which is a finding and not a reason to drop the operator.
    gap = np.zeros(n)
    np.divide(ce_expectation - ce_of_mean, ce_expectation, out=gap,
              where=covered & (ce_expectation > 0.0))
    gap_land = gap[covered]

    # What the unknown reference height is worth in the field itself, reported
    # as a bracket because it cannot be verified before a climatology exists.
    z0_bracket = [land_mean(effective_z0(k_oro, reference_height_m(t, config)))
                  for t in CE_BRACKET_K]

    output = args.output or (INPUTS / resolution.lower()
                             / f"orogen_{resolution}_surf_{ROUGHNESS_CODE:04d}.sra")
    output.parent.mkdir(parents=True, exist_ok=True)
    write_sra(output, ROUGHNESS_CODE, field)

    z_ref = {t: reference_height_m(t, config) for t in CE_BRACKET_K}

    def ce(z0v, t_air):
        """Neutral bulk exchange coefficient at the lowest model level."""
        return KARMAN ** 2 / np.log(z_ref[t_air] / np.maximum(z0v, 1e-6)) ** 2

    def bracketed(z0v):
        lo, hi = (float(ce(z0v, t)) for t in CE_BRACKET_K)
        return [round(min(lo, hi), 6), round(max(lo, hi), 6)]

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
            # SPAT-8 must hold this fixed across the ladder. It is SOLVED per
            # grid, so a difference between two rungs is part terrain and part
            # calibration: measured on the 10M mesh it moves by a factor 1.98
            # from T21 to T170 while the subgrid relief it multiplies falls by
            # 5.96. `--orographic-coefficient` supplies one value for a whole
            # ladder comparison; the default still solves, so a single-rung
            # build is unchanged.
            "coefficient_was_solved": args.orographic_coefficient is None,
        },
        "surface_aggregation": {
            "operator": "gridding.cell_expectation, NONLINEAR, over the cell's "
                        "land population: ce evaluated per mesh region and then "
                        "area-averaged, inverted to the length that reproduces "
                        "that mean. lib/gridding.py owns the reduction "
                        "operators and this file calls one rather than "
                        "reimplementing it.",
            "law": "ce = karman^2 / ln(z_ref / z0)^2",
            "population": "surface_class == land",
            "reference_height_anchor_m": round(z_ref_anchor, 2),
            "land_mean_m_over_reference_air_bracket": [round(v, 5)
                                                       for v in z0_bracket],
            # What taking `ce` of the mean length instead would have cost,
            # per cell, at the anchor height. The sign is the law's curvature:
            # `d2ce/dz0^2` carries the factor `3 - ln(z_ref/z0)`, so `ce` is
            # CONCAVE in `z0` wherever the cell is more than about twenty
            # reference heights rougher than smooth, which is nearly all land
            # here, and the expectation therefore sits BELOW `ce` of the mean.
            # A cell that is mostly playa with a canopy minority is the extreme
            # of it, and that is the case this field exists for.
            "jensen_gap_in_ce_relative": {
                "min": float(gap_land.min()),
                "median": float(np.median(gap_land)),
                "max": float(gap_land.max()),
                "land_area_weighted_mean": land_mean(gap.reshape(nlat, nlon)),
            },
            "transfer": transfer_ledger(cells, n, area, sel),
        },
        "land_mean_m": round(land_mean(z0), 5),
        "land_min_m": round(float(z0[land_cells].min()), 6),
        "land_max_m": round(float(z0[land_cells].max()), 4),
        "exchange_coefficient": {
            "note": "ce = k^2 / ln(z_ref / z0)^2, the quantity roughness acts "
                    "through. z_ref is lib/lapse.py:reference_height_m at this "
                    "planet's gravity; each value is [low, high] over "
                    "reference_air_k, across which z_ref is linear.",
            "karman": KARMAN,
            "reference_air_k": list(CE_BRACKET_K),
            "reference_height_m": [round(z_ref[t], 2) for t in CE_BRACKET_K],
            "uniform_default": bracketed(EXOPLASIM_DZ0LAND_M),
            # ce rises with z0, so the roughest cell carries the largest ce.
            "this_field_land_max": bracketed(z0[land_cells].max()),
            "this_field_land_min": bracketed(z0[land_cells].min()),
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
    lo_t, hi_t = CE_BRACKET_K
    print(f"  exchange coefficient {ce(z0[land_cells].min(), lo_t):.5f} "
          f"(flattest, barest) to {ce(z0[land_cells].max(), lo_t):.5f} "
          f"(roughest), against a uniform "
          f"{ce(EXOPLASIM_DZ0LAND_M, lo_t):.5f}, all at {lo_t} K")
    print(f"wrote {output}\n      {report_path}")


if __name__ == "__main__":
    main()
