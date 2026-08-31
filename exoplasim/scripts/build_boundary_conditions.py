#!/usr/bin/env python3
"""Write ExoPlaSim land support and topography from the World Orogen mesh.

Integrates from the mesh rather than remapping the equirectangular map PNGs.
That older path was obsolete twice over: the fork emits Gaussian grids directly
off the mesh, so a remap is a lossy intermediate, and the PNG land mask uses the
`elevation > 0` convention, which floods every dry closed-basin floor.

Two decisions worth stating, because they are where this differs from just
reading `planet.nc`:

**Land comes from `surface_class`, not `land_mask`.** The two disagree by 1.9% of
the planet, all of it dry basin floor below sea level, which is the terrain the
fork exists to preserve. `land_mask` would put it under water.

**The threshold is a rounding and its cost is reported, not hidden in a net.**
`coastline_ledger` states the two one-signed halves separately: the land a
sub-threshold cell drops into the slab ocean, and the ocean and inland water a
supra-threshold cell promotes to dry land. `surface_class` has three classes,
not two, and the threshold is applied to the LAND share alone, so a cell that
holds no ocean at all still goes to the slab ocean when its lakes take more than
half of it. See `coastline_ledger`'s own docstring for which reduction operator
produced which field.

**All fields are integrated from the native mesh, not sampled from the gridded
export.** The export resamples categorical fields, `surface_class` among them, by
taking the value of the region containing the cell centre. Against the 10M
fine-support reference a global T42 cell holds about 1,221 regions on average,
so a point sample throws away the coastline. Here the
land fraction of each cell is the area-weighted fraction of its regions that are
land. Code 1720 retains that fraction; code 172 thresholds it at
`model.geography_land_threshold` for ownership and legacy output. Topography is
the land-area-weighted mean elevation over the land regions only, so ocean
depths never drag a coastal cell's elevation down.

Elevation is written as geopotential, which is what SRA code 129 expects, using
this planet's gravity. Dry basin floors keep their negative elevation: a floor
562 m below sea level is 562 m below sea level, and flattening it to zero would
undo the preservation the whole pipeline is built around.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path

import numpy as np
import yaml

from _paths import CONFIG, INPUTS, PROJECT_ROOT
from sra import write_sra
from builds import resolution_of, grid_export, mesh_export
from gridding import cell_fraction, cell_sum, land_weighted, region_cells, transfer_ledger
from provenance import config_stamp
from orogen import Export, INLAND_WATER, LAND, OCEAN

LAND_MASK_CODE = 172
# The binary mask remains the model's ownership/topology field.  This companion
# carries the native-mesh subaerial AREA share for the two-surface-tile
# atmosphere boundary selected by SPAT-5.
LAND_FRACTION_CODE = 1720
TOPOGRAPHY_CODE = 129

TRANSPORT_LOOP = PROJECT_ROOT / "ocean" / "config" / "transport_loop.yaml"
PLANETARY_AREA_CRITERION = "sea_ice_area"


def _planetary_area_bar(path: Path = TRANSPORT_LOOP) -> dict:
    """The one bar this project declares on a planetary AREA FRACTION.

    Read from `ocean/config/transport_loop.yaml` rather than written down here,
    because a threshold copied into a script is a declaration that cannot
    propagate: nothing re-derives it and nothing objects when the two disagree.

    It is the ocean loop's `sea_ice_area` criterion, and it is the right
    instrument for a mask error because the area it was set for has the same
    climate role -- an albedo and insulation contrast switched by a threshold on
    a continuous field -- as a wet/dry partition. `notes/audits/ocean-support-
    nonlinear-reductions.md` states that reasoning and uses the same bar, so the
    ocean side of the coastline cut and the ocean support's own area errors are
    quoted in one unit.
    """
    loop = yaml.safe_load(path.read_text(encoding="utf-8"))
    for criteria in _walk_criteria(loop):
        for row in criteria:
            if isinstance(row, dict) and row.get("id") == PLANETARY_AREA_CRITERION:
                if row.get("units") != "fraction_of_planet":
                    raise ValueError(
                        f"{PLANETARY_AREA_CRITERION} is declared in "
                        f"{row.get('units')!r} and not fraction_of_planet; the "
                        "ocean-side area cost has no bar in its own units")
                return {"id": row["id"], "threshold": float(row["threshold"]),
                        "units": row["units"],
                        "source": f"{path.name}:{row['id']}",
                        "statistic": row.get("statistic"),
                        "derivation": row.get("derivation")}
    raise ValueError(f"{path} declares no {PLANETARY_AREA_CRITERION} criterion")


def _walk_criteria(node):
    """Every list named `criteria` anywhere in the loop declaration."""
    if isinstance(node, dict):
        for key, value in node.items():
            if key == "criteria" and isinstance(value, list):
                yield value
            else:
                yield from _walk_criteria(value)
    elif isinstance(node, list):
        for item in node:
            yield from _walk_criteria(item)


def coastline_ledger(mesh: Export, grid_dir: Path, mask: np.ndarray,
                     elev_m: np.ndarray, threshold_value: float) -> dict:
    """What thresholding a fractional coastline into a binary mask costs.

    The model's mask is binary and cannot be otherwise: `oceanmod.f90:256`
    hard-binarises `yls` at 0.5, so there is no tiling and no partial-water
    cell whatever this builder writes. What the threshold does is therefore not
    a choice about the model, it is a rounding, and this is its size.

    THE ROUNDING IS TWO ONE-SIGNED ERRORS THAT PARTLY CANCEL. Land in a cell
    below the threshold is dropped into the slab ocean; ocean and inland water
    in a cell above it are promoted to dry land. The report already carried the
    NET of the two, as `land_fraction_gauss_weighted` against
    `mesh_land_fraction`, and a net says nothing about either half: they can be
    several times the difference between them.

    THE THIRD CLASS IS WHY THIS IS NOT SYMMETRIC. `surface_class` partitions
    the mesh into ocean, subaerial land and inland water, and the threshold is
    applied to the LAND share alone. So inland water competes with land for its
    own cell's classification, and a cell that is entirely terrestrial -- no
    ocean in it at all -- goes to the slab ocean whenever its lakes take more
    than half of it.

    Which operator produced which field, by `lib/gridding.py`'s semantics:
    the three class shares are CATEGORICAL, `cell_fraction` over an area
    population, and they partition the cell so they sum to one wherever the
    mesh covers it; the areas and the land volume are EXTENSIVE, `cell_sum`;
    the elevation the builder writes is INTENSIVE, `cell_mean` over the land
    population, which is what `land_weighted` calls. A fractional cover is an
    AREA share and is never a share of the region COUNT.

    The land volume is the extensive closure: area times elevation about the
    datum, summed over the mesh's land, against the same sum over what the
    model receives. It is the quantity a mask change moves that an area alone
    does not, because the terrain the fork exists to preserve sits below the
    datum and enters it with the opposite sign.

    `candidate_rules` prices the ALTERNATIVES on the same grid and in the same
    terms. A threshold is a decision, and a decision taken against one rule's
    cost is not a comparison; the three that were proposed answer different
    questions -- conserve land AREA, exempt below-datum land from the threshold
    outright, or leave 0.5 -- and each is reported whether or not it is
    adopted. Nothing in that block changes what this builder writes.
    `notes/audits/coastline-threshold-cost.md` reads them.
    """
    cell, nlat, nlon = region_cells(mesh, grid_dir)
    ncell = nlat * nlon
    area = mesh.cell_area.astype(np.float64)
    sc = mesh.surface_class
    flat = mask.reshape(-1) > 0

    shares = {}
    for name, code in (("ocean", OCEAN), ("land", LAND),
                       ("inland_water", INLAND_WATER)):
        shares[name], covered = cell_fraction(cell, ncell, area, sc == code)
    total = cell_sum(cell, ncell, area)
    partition = shares["ocean"] + shares["land"] + shares["inland_water"]

    def area_of(selector, cells) -> float:
        return float(cell_sum(cell, ncell, area, selector)[cells].sum())

    is_land = sc == LAND
    is_ocean = sc == OCEAN
    mesh_land = float(area[is_land].sum())
    mesh_ocean = float(area[is_ocean].sum())
    model_land = float(total[flat].sum())
    dropped = area_of(is_land, ~flat)
    promoted_ocean = area_of(is_ocean, flat)
    promoted_inland = area_of(sc == INLAND_WATER, flat)

    # THE DENOMINATOR FOR THE OCEAN SIDE IS THE PLANET AND NOT THE MESH'S OWN
    # SUM. The bar below is a fraction of planetary area, so the area it is a
    # fraction of is the sphere the manifest declares, 4*pi*r^2, and not the sum
    # of the mesh's polygon areas, which exceeds it by about 7e-4 of the planet
    # -- the same order as the bar. Reporting both puts that difference where a
    # reader can see it rather than inside a ratio.
    planet_area = 4.0 * math.pi * mesh.radius_km ** 2
    bar = _planetary_area_bar()

    # The land the threshold drops, in the terrain the fork exists to preserve.
    below = is_land & (elev_m < 0.0)
    dropped_below = area_of(below, ~flat)
    mesh_below = float(area[below].sum())

    # Extensive closure. `elev_m` is metres and `area` is the mesh's own unit,
    # so the product is a volume in that unit times a metre; only the ratio of
    # the two routes is reported, which is what closure means here.
    mesh_volume = float((area[is_land] * elev_m[is_land]).sum())
    per_cell_elev = np.zeros(ncell)
    land_area_cell = cell_sum(cell, ncell, area, is_land)
    np.divide(cell_sum(cell, ncell, area * elev_m, is_land), land_area_cell,
              out=per_cell_elev, where=land_area_cell > 0)
    model_volume = float((total[flat] * per_cell_elev[flat]).sum())

    partial = covered & (shares["land"] > 0.0) & (shares["land"] < 1.0)

    # -- what each CANDIDATE rule would cost, on this grid ------------------
    #
    # The threshold is a decision, and a decision cannot be taken against one
    # rule's cost alone. Three rules were proposed and they answer different
    # questions: move the threshold so land AREA is conserved; exempt cells
    # holding below-datum land from the threshold, which is the only form that
    # keeps `source/README.md`'s first rule; or leave 0.5. Each is priced here,
    # in the same terms, so the comparison is read off the artifact.
    #
    # `land_share` is the quantity the threshold is applied to and `below_share`
    # is the below-datum land as a share of the cell's WHOLE area, not of its
    # land, because what an exemption has to be keyed on is how much of the cell
    # the preserved terrain actually is.
    land_share = np.zeros(ncell)
    np.divide(land_area_cell, total, out=land_share, where=covered)
    below_area_cell = cell_sum(cell, ncell, area, below)
    below_share = np.zeros(ncell)
    np.divide(below_area_cell, total, out=below_share, where=covered)

    def cost(m) -> dict:
        m = np.asarray(m, dtype=bool)
        model = float(total[m].sum())
        model_wet = float(total[covered & ~m].sum())
        ocean_gap = abs(model_wet - mesh_ocean)
        return {
            "land_cells": int(m.sum()),
            "model_land_area_relative": model / mesh_land - 1.0,
            "below_datum_land_dropped": (float(below_area_cell[~m].sum()) / mesh_below
                                         if mesh_below > 0 else 0.0),
            "land_dropped_of_mesh_land": float(
                cell_sum(cell, ncell, area, is_land)[~m].sum()) / mesh_land,
            "water_promoted_of_model_land": float(
                (cell_sum(cell, ncell, area, sc == OCEAN)[m]
                 + cell_sum(cell, ncell, area, sc == INLAND_WATER)[m]).sum()) / model,
            "land_volume_closure_relative": (
                float((total[m] * per_cell_elev[m]).sum()) - mesh_volume)
                / abs(mesh_volume),
            # THE OCEAN SIDE OF THE SAME CUT. Every row above is a share of the
            # LAND, so a rule that looks cheap against a land denominator says
            # nothing about the sea it moves; and the land side quotes no bar at
            # all, which is what stops a reader re-taking the decision. These are
            # the same three flows against the ocean's own denominators plus the
            # one declared bar on a planetary area fraction.
            "model_ocean_area_relative": model_wet / mesh_ocean - 1.0,
            "ocean_promoted_of_mesh_ocean": float(
                cell_sum(cell, ncell, area, is_ocean)[m].sum()) / mesh_ocean,
            "land_added_of_model_ocean": (
                float(cell_sum(cell, ncell, area, is_land)[covered & ~m].sum())
                / model_wet if model_wet > 0 else 0.0),
            "ocean_area_gap_km2": ocean_gap,
            "ocean_area_gap_of_planet": ocean_gap / planet_area,
            "ocean_area_gap_over_bar": ocean_gap / planet_area / bar["threshold"],
        }

    # The threshold that conserves land AREA. Bisected rather than solved: the
    # mask is a step function of the threshold, so there is no derivative.
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if float(total[land_share >= mid].sum()) > mesh_land:
            lo = mid
        else:
            hi = mid
    area_conserving = 0.5 * (lo + hi)

    candidates = {
        "why": ("a threshold is a DECISION and cannot be taken against one "
                "rule's cost alone; these are the alternatives priced in the "
                "same terms. Nothing here changes what this builder writes."),
        "as_configured": cost(mask.reshape(-1) > 0),
        "area_conserving": {"threshold": area_conserving,
                            **cost(land_share >= area_conserving)},
    }
    # The exemption family: a sub-threshold cell is kept as land when the
    # below-datum land in it reaches this share of the cell. `any` is the
    # outright exemption, which is the form that recovers ALL of the preserved
    # terrain and is priced here so what it costs to do so is visible.
    for share, name in ((1.0e-12, "any"), (0.10, "0.10"), (0.25, "0.25")):
        candidates[f"exempt_below_datum_share_{name}"] = cost(
            (land_share >= threshold_value) | (below_share >= share))

    return {
        "candidate_rules": candidates,
        "why": ("the model's mask is binary at oceanmod.f90:256 whatever this "
                "builder writes; this is the size of the rounding, in both "
                "signs, not a proposal to change it"),
        "operators": {
            "class_shares": "gridding.cell_fraction, CATEGORICAL, area share",
            "areas_and_volume": "gridding.cell_sum, EXTENSIVE",
            "elevation": "gridding.cell_mean over the land population, "
                         "INTENSIVE, via gridding.land_weighted",
        },
        "class_partition_max_residual": float(
            np.abs(partition[covered] - 1.0).max()) if covered.any() else 0.0,
        "partial_land_cells": int(partial.sum()),
        "partial_land_cell_area_fraction": float(
            total[partial].sum() / total[covered].sum()) if covered.any() else 0.0,
        "mesh_land_area": mesh_land,
        "model_land_area": model_land,
        "net_land_area_relative": (model_land - mesh_land) / mesh_land,
        # THE VERDICT'S OTHER HALF. `notes/audits/ocean-support-nonlinear-
        # reductions.md` section 2 measures this cut from the wet side and finds
        # it 3.3 to 5.0 times the bar across the ladder, on a quantity that
        # barely falls along it because a coastline does not get shorter. The
        # numbers were reconstructible from the two one-signed halves above and
        # nothing carried them, so the decision record priced one side of a
        # boundary and not the other.
        "ocean_side": {
            "why": ("the land rows are shares of the land and quote no bar; the "
                    "same cut against the ocean's own denominators, against the "
                    "one bar this project declares on a planetary area fraction"),
            "bar": bar,
            "planet_area": planet_area,
            "planet_area_source": "4*pi*r^2 from the mesh manifest's radiusKm",
            "mesh_ocean_area": mesh_ocean,
            "model_ocean_area": float(total[covered & ~flat].sum()),
            "net_ocean_area_relative": (
                float(total[covered & ~flat].sum()) - mesh_ocean) / mesh_ocean,
            "ocean_area_gap_km2": abs(float(total[covered & ~flat].sum()) - mesh_ocean),
            "ocean_area_gap_of_planet": abs(
                float(total[covered & ~flat].sum()) - mesh_ocean) / planet_area,
            "ocean_area_gap_over_bar": abs(
                float(total[covered & ~flat].sum()) - mesh_ocean)
                / planet_area / bar["threshold"],
            "goldstein_is_not_this_mask": (
                "the ocean audit's 12.3 times the bar on the 36 x 36 GOLDSTEIN "
                "candidate is that grid's own binarisation and is OCN-11's; what "
                "is priced here is the ladder ExoPlaSim binarises"),
        },
        "land_dropped_to_ocean": {
            "area": dropped,
            "of_mesh_land": dropped / mesh_land,
            "below_datum_area": dropped_below,
            "of_mesh_land_below_datum": (dropped_below / mesh_below
                                         if mesh_below > 0 else 0.0),
        },
        "water_promoted_to_land": {
            "ocean_area": promoted_ocean,
            "inland_water_area": promoted_inland,
            "of_model_land": (promoted_ocean + promoted_inland) / model_land,
        },
        "wholly_terrestrial_cells_sent_to_ocean": int(
            (~flat & covered & (shares["ocean"] <= 0.0)
             & (shares["land"] > 0.0)).sum()),
        "land_volume_closure_relative": ((model_volume - mesh_volume)
                                         / abs(mesh_volume)),
        "transfer": transfer_ledger(cell, ncell, area),
    }


def build(mesh: Export, grid_dir: Path, threshold: float, gravity: float):
    """Land mask and geopotential on `grid_dir`'s grid, integrated from the mesh.

    Elevation is metres, so the land-weighted mean is the mean elevation of the
    land in each cell. Ocean depths never enter it.
    """
    elev_m = mesh.elevation_km.astype(np.float64) * 1000.0
    fraction, mean_elev, empty = land_weighted(mesh, grid_dir, elev_m)
    mask = (fraction >= threshold).astype(np.float64)
    mean_elev = np.where(mask > 0, mean_elev, 0.0)
    return {
        "land_mask": mask,
        "cells_below_mesh_resolution": empty,
        "land_fraction": fraction,
        "geopotential": mean_elev * gravity,
        "elevation_m": mean_elev,
        "coastline_ledger": coastline_ledger(mesh, grid_dir, mask, elev_m,
                                            threshold),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--mesh", type=Path, default=None,
                    help="export carrying raw/; defaults to the configured build")
    ap.add_argument("--grid", type=Path, default=None,
                    help="export whose grid to target; defaults to the config resolution")
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    model, planet = config["model"], config["planet"]
    nlat, nlon = int(model["latitudes"]), int(model["longitudes"])
    threshold = float(model["geography_land_threshold"])
    gravity = float(planet["gravity_m_s2"])

    # Deliberately from the grid, not from config: see builds.resolution_of.
    grid_dir = args.grid or grid_export(config)
    # From the grid, not from config: the two differ exactly when someone
    # builds for another resolution, which is when the filename matters.
    resolution = resolution_of(grid_dir)
    output = args.output or (INPUTS / resolution.lower())
    ex = Export(args.mesh or mesh_export(config))
    out = build(ex, grid_dir, threshold, gravity)
    # Against the grid's own dimensions, not the config's. Config's latitudes and
    # longitudes describe the resolution a *run* uses; this script writes inputs
    # for whichever grid it was given, and refusing a valid non-default grid on
    # that basis is the check misfiring rather than catching anything. The grid
    # export states its own shape, so compare against that.
    import json as _json
    grid_shape = None
    try:
        gm = _json.loads((grid_dir / "manifest.json").read_text(encoding="utf-8"))
        grid_shape = (int(gm["grid"]["height"]), int(gm["grid"]["width"]))
    except Exception:
        pass
    if grid_shape is not None and out["land_mask"].shape != grid_shape:
        raise RuntimeError(
            f"integrated grid is {out['land_mask'].shape}, but "
            f"{grid_dir.name} declares {grid_shape}"
        )
    if grid_dir == grid_export(config) and out["land_mask"].shape != (nlat, nlon):
        raise RuntimeError(
            f"the configured grid is {out['land_mask'].shape}, config asks for "
            f"{(nlat, nlon)}; these must agree for the default grid"
        )

    output.mkdir(parents=True, exist_ok=True)
    write_sra(output / f"orogen_{resolution}_surf_{LAND_MASK_CODE:04d}.sra",
              LAND_MASK_CODE, out["land_mask"])
    write_sra(output / f"orogen_{resolution}_surf_{LAND_FRACTION_CODE:04d}.sra",
              LAND_FRACTION_CODE, out["land_fraction"])
    write_sra(output / f"orogen_{resolution}_surf_{TOPOGRAPHY_CODE:04d}.sra",
              TOPOGRAPHY_CODE, out["geopotential"])

    gw = np.fromfile(grid_dir / "grid" / "gauss_weights.bin", dtype="float64")
    m, e = out["land_mask"], out["elevation_m"]
    land_cells = m > 0
    report = {
        # Top level as well as inside `mesh`, because the consistency checker
        # looks here and a report that does not state its terrain cannot be
        # checked against the build it claims to describe.
        "terrain_hash": ex.terrain_hash,
        "generated": datetime.now(timezone.utc).isoformat(),
        "mesh": ex.provenance(),
        "grid": str(grid_dir),
        "resolution": resolution,
        "method": ("Area-weighted land fraction per cell from the native mesh, "
                   "retained on code 1720 and thresholded only for code 172; "
                   "topography is the land-area-weighted mean elevation over "
                   "land regions only."),
        "land_definition": "surface_class == land (subaerial)",
        "land_threshold": threshold,
        "gravity_m_s2": gravity,
        "land_fraction_gauss_weighted": float(
            (m.mean(axis=1) * gw).sum() / gw.sum()),
        "mesh_land_fraction": float(
            ex.cell_area[ex.surface_class == LAND].sum() / ex.cell_area.sum()),
        "land_cells": int(land_cells.sum()),
        "land_fraction": {
            "surface_code": LAND_FRACTION_CODE,
            "minimum": float(out["land_fraction"].min()),
            "maximum": float(out["land_fraction"].max()),
            "partial_cells": int(((out["land_fraction"] > 0.0)
                                  & (out["land_fraction"] < 1.0)).sum()),
            "gauss_weighted": float(
                (out["land_fraction"].mean(axis=1) * gw).sum() / gw.sum()),
        },
        "cells_below_mesh_resolution": out["cells_below_mesh_resolution"],
        # What the binary threshold costs, in both signs. SPAT-5.
        "coastline_ledger": out["coastline_ledger"],
        "elevation_m": {
            "min": float(e[land_cells].min()),
            "max": float(e[land_cells].max()),
            "mean": float(e[land_cells].mean()),
            "cells_below_sea_level": int((e[land_cells] < 0).sum()),
        },
        "codes": [TOPOGRAPHY_CODE, LAND_MASK_CODE, LAND_FRACTION_CODE],
    }
    # Provenance stamp; lib/provenance.py owns the shape and the inert set.
    report.update(config_stamp(config, "exoplasim/scripts/build_boundary_conditions.py"))
    (output / "boundary_conditions_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
